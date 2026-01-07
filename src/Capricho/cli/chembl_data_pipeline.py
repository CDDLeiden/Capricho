from inspect import signature
from pathlib import Path
from typing import Literal, Optional, Union

import pandas as pd
from chemFilters.chem.standardizers import ChemStandardizer, InchiHandling
from job_tqdflex import ParallelApplier

from ..chembl.data_flag_functions import (
    flag_censored_activity_comment,
    flag_insufficient_assay_overlap,
    flag_inter_document_duplication,
    flag_max_assay_size,
    flag_min_assay_size,
    flag_missing_canonical_smiles,
    flag_missing_document_date,
    flag_missing_standard_smiles,
    flag_patent_source,
    flag_salt_or_solvent_removal,
    flag_strict_mutant_assays,
    flag_to_remove_mixture_compounds,
    flag_undefined_stereochemistry,
    flag_unit_conversion,
)
from ..chembl.exceptions import BioactivitiesNotFoundError
from ..chembl.processing import get_bioactivities_workflow
from ..chembl.unit_conversions import (
    convert_dose_units,
    convert_mass_concentration_units,
    convert_molar_concentration_units,
    convert_permeability_units,
    convert_time_units,
)
from ..core.default_fields import (
    ASSAY_ID,
    DATA_DROPPING_COMMENT,
    DATA_PROCESSING_COMMENT,
    MOLECULE_ID,
    TARGET_ID,
)
from ..core.fp_utils import calculate_mixed_FPs
from ..core.pandas_helper import save_dataframe
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import process_repeat_mols, repeated_indices_from_array_series
from ..core.stereo import find_undefined_stereocenters
from ..logger import logger

# when aggregated, some `activity_id` values will be strings and sorting won't work properly
AGGREGATE_SAVE_SORTED_BY = ["target_chembl_id", "assay_chembl_id"]
# after the workflow, `activity_id` is an integer, so we can sort by it to ensure consistent
# ordering on the aggregated datapoints -> assay1|assay2|...|assayN,activity_id1|...|activity_idN
WORKFLOW_SAVE_SORTED_BY = [*AGGREGATE_SAVE_SORTED_BY, "activity_id"]


def _warn_info_post_aggregation_repeats(
    df: pd.DataFrame,
    extra_id_cols: list[str],
    aggregate_mutants: bool = False,
    value_col: str = "pchembl_value",
    _limit: int = 15,  # limit in the string length for the warning/info logging
) -> None:
    def _truncate_dataframe(df: pd.DataFrame, limit: int) -> pd.DataFrame:
        """Truncate DataFrame values to a specified length."""
        if pd.__version__ > "2.1.0":  # applymap got deprecated in 2.1.0
            return df.map(lambda x: str(x)[:limit] + "..." if len(str(x)) > limit else str(x))
        else:
            return df.applymap(lambda x: str(x)[:limit] + "..." if len(str(x)) > limit else str(x))

    if aggregate_mutants:
        col_subset_dupli_warning = ["connectivity", "target_chembl_id", *extra_id_cols]
    else:
        col_subset_dupli_warning = ["connectivity", "mutation", "target_chembl_id", *extra_id_cols]

    value_mean_col = f"{value_col}_mean"
    logging_subset = [  # subset of columns to be displayed on the warning/info logging
        *col_subset_dupli_warning,
        "molecule_chembl_id",
        "assay_chembl_id",
        value_mean_col,
    ]

    # Based on the ID columns, we shouldn't have any duplicates. This warning is a safeguard
    duplics_for_warning = df.duplicated(subset=col_subset_dupli_warning)
    if duplics_for_warning.any():
        dupli_subset = (
            df[duplics_for_warning]
            .loc[:, logging_subset]
            .sort_values(
                by=["target_chembl_id", "connectivity", value_mean_col],
                ascending=[True, True, False],
            )
        )
        truncated_df = _truncate_dataframe(dupli_subset, _limit)
        logger.warning(
            f"There two or more compounds matching the ID columns {col_subset_dupli_warning} "
            "This is not intentional, please further inspect the collected dataset. Here's a sample "
            "of the repeated entries:\n"
            f"{truncated_df.head(10).to_string(index=False)}"
        )

    # Additional safeguard to ensure proper handling of the output by the user prior to modeling
    target_cpd_col_subset = ["target_chembl_id", "connectivity"]
    duplics_for_info = df.duplicated(subset=target_cpd_col_subset)
    if duplics_for_info.any():
        dupli_subset = (
            df[duplics_for_info]
            .loc[:, logging_subset]
            .sort_values(by=target_cpd_col_subset + [value_mean_col], ascending=[True, True, False])
        )
        truncated_df = _truncate_dataframe(dupli_subset, _limit)
        logger.info(
            "There are two or more repeated compound-target readouts (based on `connectivity` & `target_chembl_id`) "
            "without considering other ID columns. This is a result of your aggregation criteria. Make "
            "sure to differ these data points in your modeling pipeline by including information of your other id_columns, "
            "or resolve these compound-target repeats prior to modeling. Here's a sample of the repeated entries:\n"
            f"{truncated_df.head(10).to_string(index=False)}"
        )

    return


def get_standardize_and_clean_workflow(
    molecule_ids: Optional[list[str]] = None,
    target_ids: Optional[list[str]] = None,
    assay_ids: Optional[list[str]] = None,
    document_ids: Optional[list[str]] = None,
    chirality: bool = True,
    calculate_pchembl: bool = False,
    output_path: Optional[Union[str, Path]] = None,
    confidence_scores: list[str] = [7, 8, 9],
    bioactivity_type: Optional[list[str]] = None,
    standard_relation: list[str] = ["="],
    standard_units: Optional[list[str]] = None,
    assay_types: list[str] = ["B", "F"],
    chembl_release: Optional[int] = None,
    save_not_aggregated: bool = True,
    drop_unassigned_chiral: bool = False,
    version: Optional[Union[int, str]] = None,
    backend: Literal["downloader", "webresource"] = "downloader",
    curate_annotation_errors: bool = True,
    require_doc_date: bool = False,
    min_assay_size: Optional[int] = None,
    max_assay_size: Optional[int] = None,
    min_assay_overlap: int = 0,
    strict_mutant_removal: bool = False,
    value_col: str = "pchembl_value",
    enable_unit_conversion: bool = False,
) -> pd.DataFrame:  # Changed return type annotation to pd.DataFrame
    """Fetched the filtered data from ChEMBL based on the provided IDs, assay confidence,
    and bioactivity types. The fetched smiles are then standardized and chemical mixtures
    are removed from the dataset. Duplicate data is also removed and the remaining data
    is saved to a csv file.

    Args:
        molecule_ids: list of ChEMBL molecule IDs to filter data from
        target_ids: list of ChEMBL target IDs to filter data from
        assay_ids: list of ChEMBL assay IDs to filter data from
        document_ids: list of ChEMBL document IDs to filter data from
        calculate_pchembl: whether to calculate pchembl values when not found for assay
            results reported in nanomolar/micromolar units
        chirality: setting this to False will remove stereochemistry information from the
            SMILES on top of the standardization. Defaults to True.
        output_path: path to save the resulting csv file
        confidence_scores: list of confidence scores (assay-related) to filter data from
        bioactivity_type: list of bioactivity types (assay-related) to filter data from
        standard_relation: standard relation to filter data from. Currently only supports "="
            Defaults to "=".
        chembl_release: latest ChEMBL release to retrieve data from
        save_not_aggregated: whether to save the resulting data to the csv (output_path) before
        drop_unassigned_chiral: whether to drop data points with undefined stereocenters. Defaults to False.
        version: `backend=="downloader"` only! version of the ChEMBL database to be downloaded by
            chembl_downloader. If left as None, the latest version will be downloaded. Defaults to None.
        backend: the backend to be used for fetching the data. If downloader, the ChEMBL sql database
            is downloaded and extracted first. Defaults to "downloader".
        curate_annotation_errors: Whether to apply activity curation based on pChEMBL values diverging
            in exactly 3.0 (indicate possible annotation errors). Defaults to True.
        require_doc_date: Whether to filter out activities without a document year.
        max_assay_size: Minimum number of compounds in an assay. Assays smaller than this size will
            have their activities flagged for removal. Defaults to None (no filtering).
        max_assay_size: Maximum number of compounds in an assay. Assays exceeding this size will
            have their activities flagged for removal. Defaults to None (no filtering).
        min_assay_overlap: Minimum number of overlapping compounds between two assays for the same target
            for their activities to be considered. Defaults to 0 (no filtering).
        strict_mutant_removal: If True, assays with 'mutant', 'mutation', or 'variant' in their
            description will be flagged for removal. Defaults to False.

    Returns:
        pd.DataFrame: the filtered, standardized, and cleaned data
    """
    if output_path is not None:
        if isinstance(output_path, str):
            output_path = Path(output_path)

    # -log | log transformed values reported as Log XC50, -Log XC50, etc, might not
    # have a pchembl value, but *could* still be used.  If standard_type contains `Log`,
    # the standard_value will be transferred to pchembl_value.
    if bioactivity_type is None:
        # No filter on standard_type - fetch all
        biotypes = None
    elif calculate_pchembl:
        biotypes = []
        for act in bioactivity_type:
            biotypes.extend([f"Log {act}", f"-Log {act}", act])
    else:
        if standard_relation != ["="]:
            logger.error(
                "pchembl_values are only calculated for standard_relation='='. If you want "
                "to use censored data, please set `calculate_pchembl` to True with the flag "
                "--calculate-pchembl."
            )
        biotypes = bioactivity_type

    # get_bioactivities_workflow -> fetch with either webresource or downloader -> (minimally) process bioactivities
    # -> standardization is done here -> curate bioactivity errors (if curate_annotation_errors=True)
    full_df = get_bioactivities_workflow(
        molecule_chembl_ids=molecule_ids or None,
        target_chembl_ids=target_ids or None,
        assay_chembl_ids=assay_ids or None,
        document_chembl_ids=document_ids or None,
        confidence_scores=confidence_scores,
        assay_types=assay_types,
        standard_relation=standard_relation,
        standard_type=biotypes,
        standard_units=standard_units,
        calculate_pchembl=calculate_pchembl,
        curate_annotation_errors=curate_annotation_errors,
        require_document_date=require_doc_date,
        chembl_release=chembl_release,
        version=version,
        backend=backend,
        value_col=value_col,
    )

    # Flag activities without document dates for transparency
    # Note: if require_doc_date=True, these will be hard-filtered in process_bioactivities
    full_df = flag_missing_document_date(full_df)

    # Correct censored activity comments (inactive/inconclusive) with incorrect standard_relation='='
    full_df = flag_censored_activity_comment(full_df)

    # Flag activities from patent sources for transparency
    full_df = flag_patent_source(full_df)

    # Convert units if requested
    if enable_unit_conversion:
        logger.info("Converting units to standard formats")

        full_df = convert_permeability_units(  # Convert to 10^-6 cm/s
            full_df,
            value_col=value_col,
            unit_col="standard_units",
        )
        full_df = flag_unit_conversion(full_df)

        full_df = convert_molar_concentration_units(  # Convert to nM
            full_df,
            value_col=value_col,
            unit_col="standard_units",
        )
        full_df = flag_unit_conversion(full_df)

        full_df = convert_mass_concentration_units(  # Convert to ug/mL
            full_df,
            value_col=value_col,
            unit_col="standard_units",
        )
        full_df = flag_unit_conversion(full_df)

        full_df = convert_dose_units(  # Convert to mg/kg
            full_df,
            value_col=value_col,
            unit_col="standard_units",
        )
        full_df = flag_unit_conversion(full_df)

        full_df = convert_time_units(  # Conver to hr
            full_df,
            value_col=value_col,
            unit_col="standard_units",
        )
        full_df = flag_unit_conversion(full_df)

    # Filter out activities with standard_relation not in the user-selected values
    # This is important because flag_censored_activity_comment may change '=' to '<'
    if "standard_relation" in full_df.columns and standard_relation is not None:
        excluded_relations = ~full_df["standard_relation"].isin(standard_relation)
        num_excluded = excluded_relations.sum()
        if num_excluded > 0:
            logger.info(
                f"Filtering out {num_excluded} activities with standard_relation not in {standard_relation}. "
                "These activities will be flagged for removal and saved to the _removed_subset file."
            )
            # Flag these activities for removal
            full_df.loc[excluded_relations, DATA_DROPPING_COMMENT] = full_df.loc[
                excluded_relations, DATA_DROPPING_COMMENT
            ].fillna("") + (
                full_df.loc[excluded_relations, DATA_DROPPING_COMMENT]
                .apply(lambda x: "; " if x and str(x).strip() else "")
                .fillna("")
                + f"Standard relation not in selected values {standard_relation}"
            )

    # Filter assays by size
    if min_assay_size is not None or max_assay_size is not None:
        logger.info(
            f"Filtering assays by size: min={min_assay_size}, max={max_assay_size}. "
            "Assays with insufficient size will be flagged for removal."
        )
        full_df = full_df.pipe(flag_min_assay_size, min_assay_size=min_assay_size).pipe(
            flag_max_assay_size, max_assay_size=max_assay_size
        )

    # Filter by minimum assay overlap
    if min_assay_overlap > 0 and not full_df.empty:
        logger.info(f"Filtering assays based on minimum overlap of {min_assay_overlap} compounds.")
        full_df = flag_insufficient_assay_overlap(
            df=full_df,
            min_overlap=min_assay_overlap,
            molecule_col=MOLECULE_ID,
            assay_col=ASSAY_ID,
            target_col=TARGET_ID,
            comment_col=DATA_DROPPING_COMMENT,
        )

    # Columns to remove after standardization
    # Note: standard_value and standard_units are always preserved as multivalue columns
    cols_to_remove_post_standardization = [
        "type",
        "relation",
        "units",
        "value",
        "type",
        # "description",
    ]

    logger.debug(f"All fetched bioactivity types from ChEMBL: {full_df.standard_type.unique().tolist()}")
    logger.debug(f"Filtering for bioactivity types: {biotypes}")

    # drop rows without chemical structures
    no_smiles_mask = full_df.canonical_smiles.isna()
    if no_smiles_mask.any():
        _info = full_df[no_smiles_mask].iloc[:, :6]
        logger.info(f"Dropping rows with missing canonical smiles:\n{_info}")
        full_df = full_df.drop(index=_info.index).reset_index(drop=True)

    stdzer = ChemStandardizer(
        from_smi=True, n_jobs=8, verbose=False, isomeric=chirality, progress=True, chunk_size=1000
    )

    # Filter by bioactivity_type only if it's not None
    if bioactivity_type is not None:
        df = full_df.query("standard_type.isin(@bioactivity_type)")
    else:
        df = full_df.copy()

    df = (
        df
        # standardize the smiles & clean possible solvents & salts from the string
        .pipe(flag_missing_canonical_smiles)
        .assign(standard_smiles=lambda x: stdzer(x["canonical_smiles"]))
        .dropna(subset=["standard_smiles"])  # drop if no structure is found
        .pipe(flag_salt_or_solvent_removal)
        .assign(final_smiles=lambda x: x["standard_smiles"].apply(clean_mixtures))
        .drop(columns="standard_smiles")
        .rename(columns={"final_smiles": "standard_smiles"})
        .drop(columns=[c for c in cols_to_remove_post_standardization if c in full_df.columns])
        .reset_index(drop=True)
        .copy()
    )

    # make sure we don't have Nan, can result from merging pChEMBL-lacking calculated values
    df[DATA_PROCESSING_COMMENT] = df[DATA_PROCESSING_COMMENT].fillna("")

    # Raise error if df is empty after critical processing steps
    if df.empty:
        func_params = signature(get_standardize_and_clean_workflow).parameters
        local_vars = locals()
        # Filter out df from params to avoid large object in error message
        error_params = {
            k: v for k, v in local_vars.items() if k in func_params and k != "full_df" and k != "queried_df"
        }
        raise BioactivitiesNotFoundError(parameters=error_params)

    # find mixtures in the data
    mask = df["standard_smiles"].str.contains(".", regex=False)
    n_mixtures = mask.sum()
    if n_mixtures > 0:
        df = flag_to_remove_mixture_compounds(df)
        logger.info(f"Number of mixtures: {mask.sum()}")

    # Search for undefined stereocenters within the remaining data
    if drop_unassigned_chiral:  # here we have the problem with the "." SMILES
        # Use parallel processing for finding undefined stereocenters
        logger.debug(f"Finding undefined stereocenters in {len(df)} SMILES strings using parallel processing")
        applier = ParallelApplier(
            find_undefined_stereocenters,
            df["standard_smiles"].tolist(),
            n_jobs=8,  # Use 8 cores by default
            backend="loky",
            custom_desc="Find undefined stereocenters",
            logger=logger,
            chunk_size=200,
        )
        undefined_stereo_lists = applier()
        undefined_stereo_counts = [len(x) for x in undefined_stereo_lists]

        df = df.assign(undefined_stereocenters=undefined_stereo_counts).pipe(flag_undefined_stereochemistry)
        logger.trace(f'Unassigned stereocenters: {df["undefined_stereocenters"].unique().tolist()}')
        undefined_stereo_mask = df["undefined_stereocenters"] > 0
        if undefined_stereo_mask.any():
            logger.info(f"Flagging {undefined_stereo_mask.sum()} rows with undefined stereocenters.")
            logger.debug(df[undefined_stereo_mask].iloc[:, :5])
        if df.empty:
            logger.warning("All data points have been dropped due to undefined stereocenters!!")
            return pd.DataFrame()

    # Strict mutant removal based on assay_description
    if strict_mutant_removal:
        df = flag_strict_mutant_assays(df, strict_mutant_removal=True)
        if df.empty:
            logger.warning(
                "All data points have been dropped after strict mutant removal based on assay_description."
            )
            return pd.DataFrame()

    # for the duplication we try to find the same molecule identifiers (molID, SMILES) and
    # activity outcomes (targetID, organismID, standard_value, standard_relation), but reported
    # in different ChEMBL documents (different papers) so we can keep same-readouts reported
    # by two assays performed in the same paper !
    df = flag_inter_document_duplication(df).sort_values(WORKFLOW_SAVE_SORTED_BY).reset_index(drop=True)

    # This part needs to be removed prior to data aggregation. Here, we have either
    # inorganic compounds (SMILES removed from the salt removal step), mixtures (SMILES with "."), or
    # compounds with missing activity values, which are needed for the statistics
    missing_smiles_patt = r"Missing Standard SMILES|Missing SMILES|Mixture in SMILES"
    only_salt_entry_patt = r"^\.+$"  # if only salts are present, SMILES will be just "."

    # Build the query dynamically based on which columns exist and which value_col is used
    query_parts = [
        "data_dropping_comment.str.contains(@missing_smiles_patt, na=False, regex=True)",
        r"standard_smiles.str.contains(@only_salt_entry_patt, regex=True)",
    ]

    # Only filter by value_col if it exists in the dataframe
    if value_col in df.columns:
        query_parts.append(f"{value_col}.isna()")

    removed_subset = df.query(" | ".join(query_parts)).copy()
    if output_path is not None:
        suffixes = "".join(output_path.suffixes)
        new_name = output_path.stem.split(".")[0] + "_removed_subset" + suffixes
        save_dataframe(removed_subset, output_path.with_name(new_name))
    df = df.drop(index=removed_subset.index)

    if save_not_aggregated and output_path is not None:
        suffixes = "".join(output_path.suffixes)
        new_name = output_path.stem.split(".")[0] + "_not_aggregated" + suffixes
        save_dataframe(df, output_path.with_name(new_name))
    return df


def aggregate_data(
    df,
    chirality: bool,
    metadata_cols: list[str] = [],
    extra_id_cols: list[str] = [],
    extra_multival_cols: list[str] = [],
    aggregate_mutants: bool = False,
    output_path: Optional[Union[str, Path]] = None,
    compound_equality: Literal["mixed_fp", "connectivity", "smiles"] = "connectivity",
    value_col: str = "pchembl_value",
):
    """Aggregate the data obtained from ChEMBL by:
    1) Calculate fingerprints and use those to identify same-structure compounds;
    2) Identify identical arrays from fingerprints and aggregate the data;

    Aggregated data will contain the original data separated by a semicolon and calculate
    the mean, median, standard deviation, median absolute deviation, and value counts
    for the pchembl values.

    Args:

        df: dataframe output from `CompoundMapper.cli.workflow.fetch_standardize_and_clean_workflow`
        chirality: toggle chiral-sensitive fingerprints for identifying same molecules
        extra_id_cols: additional columns to use as identifiers for the aggregation. Passing
            `["assay_chembl_id"]` to this argument, for example, will only aggregate the data
            if the compound is the same and the assay is the same.
        extra_multival_cols: list of extra columns that you'd like to keep as aggregated
            values in the final dataframe. Caveat: these columns will be displayes as (str)
            separated by `;` in the final dataframe. Defaults to [].
        aggregate_mutants: if true, will aggregate data solely based on the target_chembl_id,
            regardless of the mutation flag in ChEMBL. Defaults to False.
        output_path: path to save the aggregated data
        compound_equality: How to identify same compounds in the dataset. If "mixed_fp", uses
            mixed fingerprints (ECFP4 + RDKitFP) to identify same compounds. If "connectivity",
            uses the first part of the InChI key (connectivity) to identify same compounds.
            If "smiles", uses standardized SMILES strings directly. Defaults to "connectivity".
        value_col: Column name containing the values to aggregate statistics on.
            Defaults to "pchembl_value". Use "standard_value" for non-pChEMBL data (e.g., % inhibition).

    Returns:
        pd.DataFrame: the aggregated data
    """
    current_extra_id_cols = list(extra_id_cols)  # mutable copy

    connectivity_writer = InchiHandling(
        convert_to="connectivity", n_jobs=4, progress=True, from_smi=True, chunk_size=None
    )

    # Track whether we pre-computed connectivity to avoid recalculating after aggregation
    precomputed_connectivity = None

    if compound_equality == "mixed_fp":
        fps = calculate_mixed_FPs(  # Fingerprints are calculated to identify same molecules in the dataset
            df["standard_smiles"].tolist(), n_jobs=8, morgan_kwargs={"useChirality": chirality}, chunk_size=50
        )
        df = df.assign(id_array=fps)
    elif compound_equality == "connectivity":
        connectivities = connectivity_writer(df["standard_smiles"].tolist())
        df = df.assign(id_array=connectivities)
        # Store connectivity before censored modification so we can reuse it after aggregation
        precomputed_connectivity = pd.Series(connectivities, index=df.index)
    elif compound_equality == "smiles":
        df = df.assign(id_array=df["standard_smiles"].values)
    else:
        raise ValueError(
            f"Invalid compound_equality value: {compound_equality}. "
            "Expected 'mixed_fp', 'connectivity', or 'smiles'."
        )

    # For censored measurements (!=), include relation and value in the compound identifier
    # so they are only aggregated if they have the same value AND relation
    has_censored = df["standard_relation"].ne("=").any()
    if has_censored:
        logger.info(
            "Detected censored measurements (standard_relation != '='). "
            f"These will only be aggregated if they have identical relation AND {value_col}."
        )
        # Round value to some decimal places to avoid floating point precision issues
        if value_col == "pchembl_value":
            rounded_value = df[value_col].round(2).astype(str)
        elif value_col == "standard_value":
            rounded_value = df[value_col].round(4).astype(str)
        # For censored measurements, append relation + value to the id_array
        censored_mask = df["standard_relation"] != "="
        df.loc[censored_mask, "id_array"] = (
            df.loc[censored_mask, "id_array"].astype(str)
            + "_"
            + df.loc[censored_mask, "standard_relation"]
            + "_"
            + rounded_value[censored_mask]
        )

    # Here we have a repeat index for compounds across all fetched data. Processing which repeats
    # get aggregated (e.g.: same target ID, same `extra_id_cols`, etc) is done in `process_repeat_mols`.
    repeats_idxs = repeated_indices_from_array_series(df["id_array"])

    # Build mapping from standard_smiles to connectivity before aggregation
    # (used to avoid recalculating connectivity after aggregation)
    if precomputed_connectivity is not None:
        smiles_to_connectivity = dict(zip(df["standard_smiles"], precomputed_connectivity))

    include_metadata = [
        "doc_type",
        "doi",
        "journal",
        "year",
        "chembl_release",
        *extra_multival_cols,
        DATA_DROPPING_COMMENT,
        DATA_PROCESSING_COMMENT,
    ]

    final_data = process_repeat_mols(
        df,
        repeats_idxs,
        solve_strat="keep",
        extra_id_cols=current_extra_id_cols,
        chirality=chirality,
        extra_multival_cols=include_metadata,
        aggregate_mutants=aggregate_mutants,
        value_col=value_col,
    )

    # Assign connectivity column - reuse precomputed values when available
    if precomputed_connectivity is not None:
        final_data = final_data.assign(connectivity=final_data["smiles"].map(smiles_to_connectivity))
    else:
        final_data = final_data.assign(connectivity=lambda x: connectivity_writer(x["smiles"].tolist()))

    # reorder the columns so that connectivity comes first and processing & dropping comes last
    xtra_cols = [DATA_PROCESSING_COMMENT, DATA_DROPPING_COMMENT]
    first_columns = ["connectivity", *current_extra_id_cols, "smiles"]
    last_columns = final_data.columns.difference(first_columns + xtra_cols).tolist() + xtra_cols
    cols = ["connectivity", *current_extra_id_cols, "smiles", *last_columns]

    final_data = final_data[cols].sort_values(AGGREGATE_SAVE_SORTED_BY).reset_index(drop=True)

    _warn_info_post_aggregation_repeats(
        final_data, extra_id_cols=extra_id_cols, aggregate_mutants=aggregate_mutants, value_col=value_col
    )

    if output_path is not None:
        save_dataframe(final_data, output_path)

    return final_data


def re_aggregate_data(
    df: pd.DataFrame,
    chirality: bool,
    extra_id_cols: list[str] = [],
    extra_multival_cols: list[str] = [],
    aggregate_mutants: bool = False,
    output_path: Optional[Union[str, Path]] = None,
    compound_equality: Literal["mixed_fp", "connectivity", "smiles"] = "connectivity",
) -> pd.DataFrame:
    """Re-aggregate the data obtained from the `aggregate_data` method after dataset
    explosion. Useful for exploring the effect of different `extra_id_cols` and other
    parameters.

    Args:
        df: dataframe output from `aggregate_data`
        chirality: toggle chiral-sensitive fingerprints for identifying same molecules
        extra_id_cols: additional columns to use as identifiers for the aggregation. Passing
            `["assay_chembl_id"]` to this argument, for example, will only aggregate the data
            if the compound is the same and the assay is the same.
        extra_multival_cols: list of extra columns that you'd like to keep as aggregated
            values in the final dataframe. Caveat: these columns will be displayes as (str)
            separated by `|` (pipe) in the final dataframe. Defaults to [].
        aggregate_mutants: if true, will aggregate data solely based on the target_chembl_id,
            regardless of the mutation flag in ChEMBL. Defaults to False.
        output_path: path to save the aggregated data
        compound_equality: How to identify same compounds in the dataset. If "mixed_fp",
            uses mixed fingerprints (ECFP4 + RDKitFP) to identify same compounds. If "connectivity",
            uses the first part of the InChI key (connectivity) to identify same compounds.
            If "smiles", uses standardized SMILES strings directly. Defaults to "connectivity".

    Returns:
        pd.DataFrame: the re-aggregated data
    """
    if "processed_smiles" in df.columns:
        df = df.rename(columns={"processed_smiles": "standard_smiles"})
    if compound_equality == "connectivity" and "connectivity" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'connectivity' column.")
    if "standard_smiles" not in df.columns:
        raise ValueError("Input DataFrame must contain a 'standard_smiles' column.")
    if "smiles" not in df.columns:
        raise ValueError(
            "This method expects the output from CompoundMapper's CLI, which includes a 'smiles' column."
        )

    if compound_equality == "mixed_fp":
        fps = calculate_mixed_FPs(
            df["standard_smiles"].tolist(), n_jobs=8, morgan_kwargs={"useChirality": chirality}
        )
        id_array = pd.Series(fps, index=df.index)
    elif compound_equality == "connectivity":
        id_array = df["connectivity"]
    elif compound_equality == "smiles":
        id_array = df["standard_smiles"]
    else:
        raise ValueError(
            f"Invalid compound_equality value: {compound_equality}. "
            "Expected 'mixed_fp', 'connectivity', or 'smiles'."
        )

    # For censored measurements (!=), include relation and pchembl_value in the compound identifier
    # so they are only aggregated if they have the same value AND relation
    if "standard_relation" in df.columns:
        has_censored = df["standard_relation"].ne("=").any()
        if has_censored:
            logger.info(
                "Detected censored measurements (standard_relation != '='). "
                "These will only be aggregated if they have identical relation AND pchembl_value."
            )
            # Round pchembl_value to 2 decimal places to avoid floating point precision issues
            rounded_pchembl = df["pchembl_value"].round(2).astype(str)
            # For censored measurements, append relation + value to the id_array
            censored_mask = df["standard_relation"] != "="
            id_array = id_array.copy()  # Create a copy to avoid modifying the original
            id_array.loc[censored_mask] = (
                id_array.loc[censored_mask].astype(str)
                + "_"
                + df.loc[censored_mask, "standard_relation"]
                + "_"
                + rounded_pchembl[censored_mask]
            )

    repeats_idxs = repeated_indices_from_array_series(id_array)

    include_metadata = [
        "doc_type",
        "doi",
        "journal",
        "year",
        "chembl_release",
        *extra_multival_cols,
        DATA_DROPPING_COMMENT,
        DATA_PROCESSING_COMMENT,
    ]

    for col in include_metadata:  # make sure columns exist
        if col not in df.columns:
            raise ValueError(
                f"Column '{col}' is required in the DataFrame but is missing. "
                "Please ensure that the DataFrame contains all necessary columns."
            )

    final_data = process_repeat_mols(  # recalculate the stats given new conditions
        df,
        repeats_idxs,
        solve_strat="keep",
        extra_id_cols=extra_id_cols,
        chirality=chirality,
        extra_multival_cols=include_metadata,
        aggregate_mutants=aggregate_mutants,
    )
    connectivity_writer = InchiHandling(
        convert_to="connectivity", n_jobs=8, progress=True, from_smi=True, chunk_size=50
    )
    final_data = final_data.assign(connectivity=lambda x: connectivity_writer(x["smiles"].tolist()))

    # Reorder columns as in the original aggregate_data function
    xtra_cols = [DATA_PROCESSING_COMMENT, DATA_DROPPING_COMMENT]
    for col in xtra_cols:
        if col not in final_data.columns:
            raise ValueError(
                f"Column '{col}' is required in the DataFrame but is missing. "
                "Please ensure that the DataFrame contains all necessary columns."
            )

    # reorder the columns so that connectivity comes first and processing & dropping comes last
    xtra_cols = [DATA_PROCESSING_COMMENT, DATA_DROPPING_COMMENT]
    first_columns = ["connectivity", *extra_id_cols, "smiles"]
    last_columns = final_data.columns.difference(first_columns + xtra_cols).tolist() + xtra_cols
    cols = ["connectivity", *extra_id_cols, "smiles", *last_columns]
    final_data = final_data[cols].sort_values(AGGREGATE_SAVE_SORTED_BY).reset_index(drop=True)

    _warn_info_post_aggregation_repeats(
        final_data, extra_id_cols=extra_id_cols, aggregate_mutants=aggregate_mutants
    )

    if output_path is not None:
        save_dataframe(final_data, output_path)

    return final_data
