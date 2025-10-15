"""Collection of functions to flag compounds based on specific criteria.

Not all the functions are annotating compounds to be removed from the dataset. Some are
used to annotate processing steps that occurred during the data processing pipeline, like:

- Salt/solvent removal (applied to the canonical SMILES since they're kept as-is by CompoundMapper)
- Calculated pChEMBL value (used when the this measure is absent in ChEMBL and calculated from nM ... etc readouts)
- Potential duplicates (when this one is found, CompoundMapper will keep only one of those in the dataset. If
  the user wants to investigate, please pass the `keep_duplicates` flag to True)

"""

from typing import Optional

import numpy as np
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


def flag_insufficient_assay_overlap(
    df: pd.DataFrame,
    min_overlap: int = 0,
    molecule_col: str = MOLECULE_ID,
    assay_col: str = ASSAY_ID,
    target_col: str = TARGET_ID,
    comment_col: str = DATA_DROPPING_COMMENT,
) -> pd.DataFrame:
    """Mark activities from assay pairs (for the same target) that don't meet the minimum
    compound overlap criterium, useful for analysis assessing the comparability of assays
    reported in ChEMBL. Depending on the target, this filter can remove a significant
    amount of activities from the dataset, but it is useful to assess the comparability of
    the assays reported in the database.

    Args:
        df: DataFrame to be processed.
        min_overlap: Minimum number of overlapping compounds required. Defaults to 0
        molecule_col: Name of the molecule identifier column. Defaults to molecule_chembl_id.
        assay_col: Name of the assay identifier column. Defaults to assay_chembl_id.
        target_col: Name of the target identifier column. Defaults to target_chembl_id.
        comment_col: Name of the column to store dropping comments. Defaults to data_dropping_comment.

    Returns:
        pd.DataFrame: DataFrame with activities from low-overlap assay pairs flagged.
    """
    if min_overlap is None:
        raise ValueError("min_overlap must be provided and cannot be None.")
        return df
    elif min_overlap == 0:
        logger.info("min_overlap set to 0. Skipping insufficient assay overlap filtering.")
        return df

    if not all(col in df.columns for col in [molecule_col, assay_col, target_col]):
        logger.warning(
            f"One or more required columns ({molecule_col}, {assay_col}, {target_col}) "
            "not found in DataFrame. Skipping insufficient assay overlap filtering."
        )
        return df

    if df.empty:
        logger.info("DataFrame is empty. Skipping insufficient assay overlap filtering.")
        return df

    if comment_col not in df.columns:  # Ensure comment_col exists
        df[comment_col] = pd.Series(dtype="object")

    assays_to_flag_globally = set()  # Keep track of assays that are part of any failing pair

    for target_id, group_df in df.groupby(target_col):
        unique_assays_in_target = group_df[assay_col].unique()
        if len(unique_assays_in_target) < 2:
            continue  # Not enough assays for this target to form a pair

        # Create a mapping of assay_id to set of molecule_ids for this target
        assay_to_molecules = {}
        for assay_id in unique_assays_in_target:
            assay_to_molecules[assay_id] = set(
                group_df[group_df[assay_col] == assay_id][molecule_col].unique()
            )

        # Iterate through unique pairs of assays for the current target
        from itertools import combinations

        for assay1_id, assay2_id in combinations(unique_assays_in_target, 2):
            molecules1 = assay_to_molecules.get(assay1_id, set())
            molecules2 = assay_to_molecules.get(assay2_id, set())

            overlap_count = len(molecules1.intersection(molecules2))

            if overlap_count < min_overlap:
                assays_to_flag_globally.add(assay1_id)
                assays_to_flag_globally.add(assay2_id)
                logger.debug(
                    f"Target {target_id}: Assays {assay1_id} and {assay2_id} have overlap {overlap_count} "
                    f"(less than {min_overlap}). They will be flagged."
                )

    if assays_to_flag_globally:
        filter_mask = df[assay_col].isin(list(assays_to_flag_globally))
        num_activities_flagged = filter_mask.sum()
        num_assays_flagged = len(assays_to_flag_globally)

        logger.info(
            f"Flagging {num_activities_flagged} activities from {num_assays_flagged} unique assays "
            f"that were part of pairs not meeting the minimum overlap of {min_overlap} compounds."
        )

        comment_text = f"Insufficient assay overlap (min_overlap={min_overlap})"

        existing_comments = df.loc[filter_mask, comment_col].fillna("")
        new_comments = existing_comments.apply(lambda x: f"{x}; {comment_text}" if x else comment_text)
        df.loc[filter_mask, comment_col] = new_comments
    else:
        logger.info(f"No assay pairs found below the minimum overlap of {min_overlap} compounds.")

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
            comment="Corrected standard_relation from '=' to '<' (censored activity_comment)",
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
