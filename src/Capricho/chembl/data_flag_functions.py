"""Collection of functions to flag compounds based on specific criteria.

Not all the functions are annotating compounds to be removed from the dataset. Some are
used to annotate processing steps that occurred during the data processing pipeline, like:

- Salt/solvent removal (applied to the canonical SMILES since they're kept as-is by CompoundMapper)
- Calculated pChEMBL value (used when the this measure is absent in ChEMBL and calculated from nM ... etc readouts)
- Potential duplicates (when this one is found, CompoundMapper will keep only one of those in the dataset. If
  the user wants to investigate, please pass the `keep_duplicates` flag to True)

"""

from typing import Optional

import pandas as pd

from ..core.default_fields import (
    ASSAY_ID,
    DATA_DROPPING_COMMENT,
    MOLECULE_ID,
    TARGET_ID,
)
from ..core.pandas_helper import add_comment, conflicting_duplicates
from ..logger import logger


### Marking readouts that will be DROPPED by CompoundMapper ###
def flag_missing_canonical_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'canonical_smiles' is missing."""
    return add_comment(
        df,
        comment="Missing SMILES",
        criteria_func=lambda x: x.isna(),  # x is the Series df["canonical_smiles"]
        target_column="canonical_smiles",
        comment_type="d",
    )


def flag_missing_standard_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'standard_smiles' is missing."""
    return add_comment(
        df,
        comment="Missing Standard SMILES",
        criteria_func=lambda x: x.isna(),
        target_column="standard_smiles",
        comment_type="d",
    )


def flag_with_data_validity_comment(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'data_validity_comment' is present (not NA)."""
    return add_comment(
        df,
        comment="Data Validity Comment Present",
        criteria_func=lambda x: ~x.isna(),
        target_column="data_validity_comment",
        comment_type="d",
    )


def flag_potential_duplicate(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'potential_duplicate' is 0."""
    return add_comment(
        df,
        comment="Potential Duplicate",
        criteria_func=lambda x: x == 1,
        target_column="potential_duplicate",
        comment_type="d",
    )


def flag_to_remove_mixture_compounds(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'mixture_compounds' is True."""
    return add_comment(
        df,
        comment="Mixture in SMILES",
        criteria_func=lambda x: x.str.contains(".", regex=False),
        target_column="standard_smiles",
        comment_type="d",
    )


def flag_undefined_stereochemistry(df: pd.DataFrame) -> pd.DataFrame:
    """Mark compounds with undefined stereochemistry based on a predefined boolean mask."""
    return add_comment(
        df,
        comment="Undefined Stereochemistry",
        criteria_func=lambda x: x > 0,
        target_column="undefined_stereocenters",
        comment_type="d",
    )


def flag_min_assay_size(df: pd.DataFrame, min_assay_size: int = 0) -> pd.DataFrame:
    """Mark assays for removal based on size lower than the specified minimum assay size."""
    assert min_assay_size >= 0, "Minimum assay size must be a non-negative integer."

    if "assay_size" not in df.columns:
        logger.error("Column 'assay_size' not found in DataFrame. Skipping minimum assay size filtering.")
        return df
    if min_assay_size == 0:
        logger.info("Minimum assay size is set to 0. Skipping filtering based on assay size.")
        return df
    else:
        return add_comment(
            df,
            comment=f"Assay size < {min_assay_size}",
            criteria_func=lambda x: x < min_assay_size,
            target_column="assay_size",
            comment_type="d",
        )


def flag_max_assay_size(df: pd.DataFrame, max_assay_size: int = 1000) -> pd.DataFrame:
    """Mark assays for removal based on size greater than the specified maximum assay size."""
    assert max_assay_size > 0, "Maximum assay size must be a positive integer."

    if "assay_size" not in df.columns:
        logger.error("Column 'assay_size' not found in DataFrame. Skipping maximum assay size filtering.")
        return df
    if max_assay_size == 0:
        logger.info("Maximum assay size is set to 0. Skipping filtering based on assay size.")
        return df
    else:
        return add_comment(
            df,
            comment=f"Assay size > {max_assay_size}",
            criteria_func=lambda x: x > max_assay_size,
            target_column="assay_size",
            comment_type="d",
        )


def flag_strict_mutant_assays(df: pd.DataFrame, strict_mutant_removal: bool = False) -> pd.DataFrame:
    """Mark assays for removal if their description contains mutant-related keywords
    and strict_mutant_removal is True.
    """
    if not strict_mutant_removal:
        logger.debug("Strict mutant removal is False. Skipping assay description-based mutant flagging.")
        return df

    if "assay_description" not in df.columns:
        logger.warning(
            "Column 'assay_description' not found in DataFrame. " "Skipping strict mutant assay flagging."
        )
        return df

    keywords = ["mutant", "mutation", "variant"]
    keyword_pattern = "|".join(keywords)

    criteria = (  # noqa: E731 - lambda function saved in a variable here just for clarity
        lambda series: series.astype(str).str.lower().str.contains(keyword_pattern, regex=True, na=False)
    )

    # Calculate the mask for logging purposes before applying add_comment
    mask_to_flag = criteria(df["assay_description"])
    num_to_flag = mask_to_flag.sum()

    if num_to_flag > 0:
        logger.info(
            f"Flagging {num_to_flag} assays for removal based on strict mutant removal criteria "
            f"(keywords: {', '.join(keywords)} in 'assay_description')."
        )
        df = add_comment(
            df,
            comment="Mutation keyword in assay description",
            criteria_func=criteria,
            target_column="assay_description",
            comment_type="d",
        )
    else:
        logger.info(
            "No assays to flag for removal based on strict mutant removal criteria in 'assay_description'."
        )
    return df


def flag_missing_document_date(df: pd.DataFrame) -> pd.DataFrame:
    """Mark activities that lack a document date (year) in the processing comment.

    This function always flags missing document dates for transparency, regardless of whether
    they will be filtered out. Activities without document dates are flagged in the
    data_processing_comment column so users can see which data points lack temporal information.

    Args:
        df: DataFrame to be processed.

    Returns:
        pd.DataFrame: DataFrame with activities lacking document dates flagged in processing comment.
    """
    if "year" not in df.columns:
        logger.debug("Column 'year' not found in DataFrame. Skipping document date flagging.")
        return df

    mask = df["year"].isna()
    num_to_flag = mask.sum()

    if num_to_flag > 0:
        logger.info(f"Flagging {num_to_flag} activities with missing document date (year) for transparency.")
        df = add_comment(
            df,
            comment="Missing document date",
            criteria_func=lambda x: x.isna(),
            target_column="year",
            comment_type="d",
        )
    else:
        logger.debug("No activities with missing document dates found.")

    return df


def flag_incompatible_units(df: pd.DataFrame) -> pd.DataFrame:
    """Mark activities with units that cannot be converted to pChEMBL in the dropping comment.

    This function flags activities with standard_units that are incompatible with pChEMBL
    calculation (i.e., not nM, µM, uM, or mM). These activities will have pchembl_value=NaN
    and are flagged for transparency.

    Args:
        df: DataFrame to be processed.

    Returns:
        pd.DataFrame: DataFrame with incompatible units flagged in dropping comment.
    """
    if "standard_units" not in df.columns:
        logger.debug("Column 'standard_units' not found in DataFrame. Skipping incompatible units flagging.")
        return df

    compatible_units = ["nM", "µM", "uM", "mM"]
    # Flag rows where standard_units is NOT in compatible_units AND NOT null
    mask = ~df["standard_units"].isin(compatible_units) & df["standard_units"].notna()
    num_to_flag = mask.sum()

    if num_to_flag > 0:
        logger.info(f"Flagging {num_to_flag} activities with incompatible units for pChEMBL calculation.")
        df = add_comment(
            df,
            comment="Incompatible units for pChEMBL calculation",
            criteria_func=lambda x: ~x.isin(compatible_units) & x.notna(),
            target_column="standard_units",
            comment_type="d",
        )
    else:
        logger.debug("No activities with incompatible units found.")

    return df


def flag_patent_source(df: pd.DataFrame) -> pd.DataFrame:
    """Mark activities that originate from patent documents.

    Patent-sourced data in ChEMBL often contains annotation errors or inconsistent
    measurements compared to peer-reviewed publications. This function flags activities
    from patents to allow users to assess data quality differences between patent and
    non-patent sources.

    Args:
        df: DataFrame to be processed.

    Returns:
        pd.DataFrame: DataFrame with patent-sourced activities flagged in dropping comment.
    """
    mask = df["doc_type"] == "PATENT"
    num_to_flag = mask.sum()

    if num_to_flag > 0:
        logger.info(f"Flagging {num_to_flag} activities from patent sources.")
        df = add_comment(
            df,
            comment="Patent source",
            criteria_func=lambda x: x == "PATENT",
            target_column="doc_type",
            comment_type="d",
        )
    else:
        logger.debug("No patent-sourced activities found.")

    return df


def flag_insufficient_assay_overlap(
    df: pd.DataFrame,
    min_overlap: int = 0,
    molecule_col: str = MOLECULE_ID,
    assay_col: str = ASSAY_ID,
    target_col: str = TARGET_ID,
    comment_col: str = DATA_DROPPING_COMMENT,
    max_assay_match: bool = False,
    assay_match_fields: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Mark activities from assay pairs (for the same target) that don't meet the minimum
    compound overlap criterium, useful for analysis assessing the comparability of assays
    reported in ChEMBL. Depending on the target, this filter can remove a significant
    amount of activities from the dataset, but it is useful to assess the comparability of
    the assays reported in the database.

    This function calculates overlap across ALL assays regardless of size flags, following
    CAPRICHO's principle of transparency. Overlap is only counted when:
    1. Compounds have DIFFERENT pChEMBL values across assays (same values indicate annotation errors)
    2. The difference is not exactly 3.0 or 6.0 log units (likely censored/inactive measurements)
    3. Assays are from DIFFERENT documents (same-document overlaps are excluded)

    When max_assay_match is True, only assay pairs with matching metadata (as defined by
    assay_match_fields) are compared for overlap. Assays that don't have any compatible
    partner with sufficient overlap will be flagged for removal.

    Args:
        df: DataFrame to be processed.
        min_overlap: Minimum number of overlapping compounds required. Defaults to 0
        molecule_col: Name of the molecule identifier column. Defaults to molecule_chembl_id.
        assay_col: Name of the assay identifier column. Defaults to assay_chembl_id.
        target_col: Name of the target identifier column. Defaults to target_chembl_id.
        comment_col: Name of the column to store dropping comments. Defaults to data_dropping_comment.
        max_assay_match: If True, only compare assays with matching metadata fields. Defaults to False.
        assay_match_fields: List of metadata fields that must match for assays to be compared.
            Only used when max_assay_match is True. If None, uses a default set of fields.

    Returns:
        pd.DataFrame: DataFrame with activities from low-overlap assay pairs flagged.
    """
    if min_overlap is None:
        raise ValueError("min_overlap must be provided and cannot be None.")
        return df
    elif min_overlap == 0:
        logger.info("min_overlap set to 0. Skipping insufficient assay overlap filtering.")
        return df

    # Check for required columns
    required_cols = [molecule_col, assay_col, target_col, "pchembl_value", "document_chembl_id"]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.warning(
            f"Required columns missing: {missing_cols}. " "Skipping insufficient assay overlap filtering."
        )
        return df

    if df.empty:
        logger.info("DataFrame is empty. Skipping insufficient assay overlap filtering.")
        return df

    if comment_col not in df.columns:  # Ensure comment_col exists
        df[comment_col] = pd.Series(dtype="object")

    # Log the overlap calculation strategy
    logger.info(
        "Calculating assay overlap with the following criteria:\n"
        "  - Only counting compounds with DIFFERENT pChEMBL values across assays\n"
        "  - Excluding differences of exactly 3.0 or 6.0 log units (likely annotation errors)\n"
        "  - Excluding overlaps within the same document"
    )

    # Validate assay_match_fields if max_assay_match is enabled
    if max_assay_match:
        if assay_match_fields is None:
            from ..core.default_fields import DEFAULT_ASSAY_MATCH_FIELDS

            assay_match_fields = DEFAULT_ASSAY_MATCH_FIELDS

        # Check which match fields are available in the DataFrame
        missing_fields = [col for col in assay_match_fields if col not in df.columns]
        if missing_fields:
            logger.warning(
                f"max_assay_match enabled but the following match fields are missing from DataFrame: {missing_fields}. "
                f"These fields will be ignored for metadata matching."
            )
            assay_match_fields = [col for col in assay_match_fields if col in df.columns]

        if not assay_match_fields:
            logger.warning(
                "max_assay_match enabled but no match fields are available in DataFrame. "
                "Falling back to standard overlap checking without metadata matching."
            )
            max_assay_match = False
        else:
            logger.info(
                f"max_assay_match enabled. Assay pairs will only be compared if they match on: {assay_match_fields}"
            )

    import numpy as np

    assays_to_flag_globally = set()  # Keep track of assays that don't have a compatible partner

    for target_id, group_df in df.groupby(target_col):
        unique_assays = group_df[assay_col].unique()

        if len(unique_assays) < 2:
            continue  # Not enough assays for this target to form a pair

        # Prepare data with necessary columns for vectorized operations
        target_data = group_df[[assay_col, molecule_col, "pchembl_value", "document_chembl_id"]].copy()

        # Add metadata columns if max_assay_match is enabled
        if max_assay_match:
            target_data = target_data.merge(
                group_df[[assay_col] + assay_match_fields].drop_duplicates(subset=[assay_col]),
                on=assay_col,
                how="left",
            )

        # Self-join to create all pairs of (assay1, assay2) for the same molecule
        # This finds all overlapping molecules between assays
        pairs = target_data.merge(target_data, on=molecule_col, suffixes=("_1", "_2"))

        # Filter to only keep assay1 < assay2 to avoid duplicates and self-comparisons
        pairs = pairs[pairs[f"{assay_col}_1"] < pairs[f"{assay_col}_2"]]

        # Apply filters
        # 1. Different documents
        pairs = pairs[pairs["document_chembl_id_1"] != pairs["document_chembl_id_2"]]

        # 2. Different pChEMBL values (excluding 3.0 and 6.0 log unit differences)
        pchembl_diff = np.abs(pairs["pchembl_value_1"] - pairs["pchembl_value_2"])
        pairs = pairs[
            (pairs["pchembl_value_1"] != pairs["pchembl_value_2"])
            & (pchembl_diff != 3.0)
            & (pchembl_diff != 6.0)
        ]

        # 3. Metadata matching if required
        if max_assay_match:
            metadata_match = pd.Series(True, index=pairs.index)
            for field in assay_match_fields:
                # Fill NaN values with a sentinel so that NaN == NaN (both missing the same metadata)
                field1 = pairs[f"{field}_1"].fillna("__MISSING__")
                field2 = pairs[f"{field}_2"].fillna("__MISSING__")
                metadata_match &= field1 == field2
            pairs = pairs[metadata_match]

        # Count overlapping compounds per assay pair
        overlap_counts = pairs.groupby([f"{assay_col}_1", f"{assay_col}_2"]).size()

        # Find assay pairs with sufficient overlap
        sufficient_pairs = overlap_counts[overlap_counts >= min_overlap]

        # Collect all assays that have at least one compatible partner
        assays_with_sufficient_overlap = set()
        for (assay1, assay2), count in sufficient_pairs.items():
            assays_with_sufficient_overlap.add(assay1)
            assays_with_sufficient_overlap.add(assay2)
            logger.trace(
                f"Target {target_id}: Assays {assay1} and {assay2} have {count} compounds "
                f"with conflicting pChEMBL values (>= {min_overlap})."
            )

        # Flag assays that don't have any compatible partner
        for assay_id in unique_assays:
            if assay_id not in assays_with_sufficient_overlap:
                assays_to_flag_globally.add(assay_id)
                logger.debug(
                    f"Target {target_id}: Assay {assay_id} has no compatible partner with >= {min_overlap} "
                    f"overlapping compounds with conflicting pChEMBL values. Will be flagged for removal."
                )

    if assays_to_flag_globally:
        num_assays_flagged = len(assays_to_flag_globally)

        if max_assay_match:
            comment_text = f"Insufficient assay overlap with metadata matching (min_overlap={min_overlap})"
            logger.info(
                f"Flagging activities from {num_assays_flagged} unique assays "
                f"that lack a metadata-compatible partner with >= {min_overlap} overlapping compounds."
            )
        else:
            comment_text = f"Insufficient assay overlap (min_overlap={min_overlap})"
            logger.info(
                f"Flagging activities from {num_assays_flagged} unique assays "
                f"that lack a partner with >= {min_overlap} overlapping compounds."
            )

        # Use add_comment helper for consistency
        df = add_comment(
            df,
            comment=comment_text,
            criteria_func=lambda x: x.isin(list(assays_to_flag_globally)),
            target_column=assay_col,
            comment_type="d",
        )
    else:
        if max_assay_match:
            logger.info(
                f"All assays have at least one metadata-compatible partner with >= {min_overlap} overlapping compounds."
            )
        else:
            logger.info(f"All assays have at least one partner with >= {min_overlap} overlapping compounds.")

    return df


def flag_censored_activity_comment(df: pd.DataFrame) -> pd.DataFrame:
    """Mark activities with activity_comment indicating censored/inactive data but standard_relation='='.

    ChEMBL contains many activities where the activity_comment field indicates the compound
    was inactive, inconclusive, or not tested, but the standard_relation is incorrectly set to '='.

    The standard_value represents a concentration (e.g., IC50 in nM), and pChEMBL = -log10(standard_value_in_M).
    When a compound is marked as "inactive" at a given concentration, it means:
    - The true IC50 (standard_value) is GREATER THAN the tested concentration
    - The true pChEMBL value is LESS THAN the reported pChEMBL

    Therefore, the standard_relation should be '<' for pChEMBL values (or '>' for standard_value).
    This function corrects the relation to '<' since we work with pChEMBL values.
    """
    if "activity_comment" not in df.columns:
        logger.debug("Column 'activity_comment' not found. Skipping censored activity comment correction.")
        return df

    if "standard_relation" not in df.columns:
        logger.warning("Column 'standard_relation' not found. Cannot correct censored activity comments.")
        return df

    # Keywords indicating censored/inactive data (case-insensitive)
    censored_keywords = [
        "not active",
        "inactive",
        "inconclusive",
        "not tested",
        "not determined",
        "nd",
        "below threshold",
        "below detection",
        "no activity",
    ]

    # Build pattern for case-insensitive matching
    pattern = "|".join([f"(?i){keyword}" for keyword in censored_keywords])

    mask = (  # Find rows with problematic activity_comment and standard_relation='='
        df["activity_comment"].notna()
        & df["activity_comment"].astype(str).str.contains(pattern, regex=True, na=False)
        & (df["standard_relation"] == "=")
    )

    num_to_correct = mask.sum()

    if num_to_correct > 0:
        logger.info(
            f"Correcting {num_to_correct} activities with censored activity_comment "
            f"but standard_relation='=' (changing to '<' for log-transformed pChEMBL values)"
        )

        # Correct standard_relation to "<";
        # inactive means actual activity is lower than reported standard_value
        df.loc[mask, "standard_relation"] = "<"

        df = add_comment(
            df,
            comment="Corrected standard_relation from = to < (censored activity_comment)",
            criteria_func=lambda x, idxs=df[mask].index: x.index.isin(idxs),
            target_column="activity_comment",
            comment_type="p",
        )
    else:
        logger.debug("No activities with censored activity_comment and standard_relation='=' found.")

    return df


### Marking readouts that were PROCESSED by CompoundMapper ###
def flag_calculated_pchembl(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'calculated_pchembl' is True."""
    return add_comment(
        df,
        comment="Calculated pChEMBL",
        comment_type="p",
    )


def flag_salt_or_solvent_removal(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows with a salt/mixture on the canonical SMILES (not modified by CompoundMapper)"""
    return add_comment(
        df,
        comment="Salt/solvent removed",
        criteria_func=lambda x: x.str.contains(".", regex=False),
        target_column="canonical_smiles",
        comment_type="p",
    )


def flag_inter_document_duplication(
    df: pd.DataFrame,
    key_subset: list[str] = [
        "molecule_chembl_id",
        "standard_smiles",
        "canonical_smiles",
        "pchembl_value",
        "standard_relation",
        "target_chembl_id",
        "mutation",
        "target_organism",
    ],
    diff_subset: Optional[list[str]] = ["document_chembl_id"],
) -> pd.DataFrame:
    """Marks rows with a potential duplication after SMILES standardization & salt removal.

    This function only flags duplicates for discrete measurements (standard_relation='=').
    Censored measurements (e.g., '<', '>') are not flagged as duplicates since the same
    bound can be independently reached in different studies without indicating true duplication.

    Args:
        df: DataFrame to be processed.
        key_subset: metadata columns used for identifying duplicates. Defaults to a
            list of columns typically used to identify a compound readout.
        diff_subset: optional metadata columns used to identify duplicates across
            different documents. If None, it identifies only based on `key_subset`.
            Defaults to a list containing 'document_chembl_id', which is used to identify
            duplicates across different documents.

    Returns:
        pd.DataFrame: DataFrame with duplicates marked in `data_processing_comment`
                      (or `data_dropping_comment` if comment_type='d' was used).
    """
    # Only check for duplicates in discrete measurements (standard_relation='=')
    if "standard_relation" not in df.columns:
        logger.warning(
            "Column 'standard_relation' not found in DataFrame. Skipping inter-document duplication flagging."
        )
        return df

    discrete_mask = df["standard_relation"] == "="
    num_discrete = discrete_mask.sum()
    num_censored = (~discrete_mask).sum()

    if num_discrete == 0:
        logger.info(
            f"No discrete measurements (standard_relation='=') found. "
            f"All {num_censored} measurements are censored. Skipping inter-document duplication flagging."
        )
        return df

    logger.debug(
        f"Checking for inter-document duplicates among {num_discrete} discrete measurements. "
        f"Ignoring {num_censored} censored measurements."
    )

    # Only check duplicates among discrete measurements
    df_discrete = df[discrete_mask]
    dupli_mask_discrete = conflicting_duplicates(df_discrete, key_subset=key_subset, diff_subset=diff_subset)
    n_duplicates = dupli_mask_discrete.sum()

    if n_duplicates > 0:
        logger.info(
            f"Flagged {n_duplicates} duplicates with same Mol identifiers across different Documents "
            f"(only among discrete measurements with standard_relation='=')."
        )
        # Map the duplicate indices back to the original dataframe
        dupli_idxs = df_discrete[dupli_mask_discrete].index
        return (
            df.assign(temp_dupli_flag=False)
            .pipe(lambda x: x.assign(temp_dupli_flag=x.index.isin(dupli_idxs)))
            .pipe(
                add_comment,
                comment="pChEMBL Duplication Across Documents",
                criteria_func=lambda x, idxs=dupli_idxs: x.index.isin(idxs),
                target_column="temp_dupli_flag",
                comment_type="p",
            )
            .drop(columns=["temp_dupli_flag"])
        )
    else:
        logger.debug("No inter-document duplicates found among discrete measurements.")
        return df


def flag_unit_conversion(df: pd.DataFrame) -> pd.DataFrame:
    """Mark rows where unit conversion was applied to standardize measurements.

    This function flags activities that had their standard_value and standard_units
    converted to a common unit by unit conversion functions (e.g., convert_permeability_units,
    convert_molar_concentration_units, etc.). The conversion_factor column (added during
    conversion) is used to identify which rows were converted.

    If an 'original_unit' column exists (added by newer conversion functions), it will be
    used to create a more informative comment showing the original -> target unit transformation.
    Both conversion_factor and original_unit columns are removed after flagging.

    This is a processing flag (comment_type='p') to document data transformations
    for transparency.

    Args:
        df: DataFrame to be processed. Must contain 'conversion_factor' column
            if unit conversion was applied. May optionally contain 'original_unit'
            for more detailed flagging.

    Returns:
        pd.DataFrame: DataFrame with converted rows flagged in data_processing_comment.
            The conversion_factor and original_unit columns are removed after flagging.
    """
    if "conversion_factor" not in df.columns:
        logger.debug("Column 'conversion_factor' not found. Skipping unit conversion flagging.")
        return df

    mask = df["conversion_factor"].notna()
    num_converted = mask.sum()

    if num_converted > 0:
        logger.info(f"Flagging {num_converted} activities with converted units.")

        # Check if we have original_unit and standard_units for detailed comment
        has_original_unit = "original_unit" in df.columns
        has_standard_units = "standard_units" in df.columns

        if not has_original_unit or not has_standard_units:
            logger.warning(
                "Unit conversion flagging requires 'original_unit' and 'standard_units' columns. "
                "Skipping flagging. This may indicate a bug in the conversion function."
            )
        else:
            # Create row-specific comments showing original -> target unit
            for idx in df[mask].index:
                original = df.loc[idx, "original_unit"]
                target = df.loc[idx, "standard_units"]
                comment = f"Unit converted to {target} from {original}"
                df = add_comment(
                    df,
                    comment=comment,
                    criteria_func=lambda x, i=idx: x.index == i,
                    target_column="conversion_factor",
                    comment_type="p",
                )
    else:
        logger.debug("No activities with converted units found.")

    # Clean up conversion_factor column
    df = df.drop(columns=["conversion_factor"])

    # Clean up original_unit column if it exists
    if "original_unit" in df.columns:
        df = df.drop(columns=["original_unit"])

    return df
