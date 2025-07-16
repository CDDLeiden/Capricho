from inspect import signature
from pathlib import Path
from typing import Literal, Optional, Union

import pandas as pd
from chemFilters.chem.standardizers import ChemStandardizer, InchiHandling

from ..chembl.data_flag_functions import (
    flag_insufficient_assay_overlap,
    flag_inter_document_duplication,
    flag_max_assay_size,
    flag_min_assay_size,
    flag_missing_canonical_smiles,
    flag_missing_standard_smiles,
    flag_salt_or_solvent_removal,
    flag_strict_mutant_assays,
    flag_to_remove_mixture_compounds,
    flag_undefined_stereochemistry,
)
from ..chembl.exceptions import BioactivitiesNotFoundError
from ..chembl.processing import get_bioactivities_workflow
from ..core.default_fields import (
    ASSAY_ID,
    DATA_DROPPING_COMMENT,
    DATA_PROCESSING_COMMENT,
    DEFAULT_ASSAY_MATCH_FIELDS,
    MOLECULE_ID,
    TARGET_ID,
)
from ..core.fp_utils import calculate_mixed_FPs
from ..core.pandas_helper import save_dataframe
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import process_repeat_mols, repeated_indices_from_array_series
from ..core.stereo import find_undefined_stereocenters
from ..logger import logger


def get_standardize_and_clean_workflow(
    molecule_ids: list[str],
    target_ids: list[str],
    assay_ids: list[str],
    document_ids: list[str],
    chirality: bool = True,
    calculate_pchembl: bool = False,
    output_path: Optional[Union[str, Path]] = None,
    confidence_scores: list[str] = [7, 8, 9],
    bioactivity_type: list[str] = ["IC50", "EC50", "AC50", "Ki", "Kd"],
    standard_relation: list[str] = ["="],  # TODO: later support data with  >, <, >=, <=...
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
                for their activities to be considered. Defaults to None (no filtering).
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
    # TODO: I noticed some negative values reported in standard_value, maybe it's because
    # the standard_type was Log? If so, we could try rescuing those values to pchembl_value
    if calculate_pchembl:
        biotypes = []
        for act in bioactivity_type:
            biotypes.extend([f"Log {act}", f"-Log {act}", act])
    else:
        biotypes = bioactivity_type

    full_df = get_bioactivities_workflow(
        molecule_chembl_ids=molecule_ids,
        target_chembl_ids=target_ids,
        assay_chembl_ids=assay_ids,
        document_chembl_ids=document_ids,
        confidence_scores=confidence_scores,
        assay_types=assay_types,
        calculate_pchembl=calculate_pchembl,
        curate_annotation_errors=curate_annotation_errors,
        require_document_date=require_doc_date,
        chembl_release=chembl_release,
        version=version,
        backend=backend,
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

    cols_to_remove_post_standardization = [
        "type",
        "relation",
        "units",
        "value",
        "standard_value",  # we'll use pchembl instead
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

    stdzer = ChemStandardizer(from_smi=True, n_jobs=8, verbose=False, isomeric=chirality)
    df = (
        full_df.query("standard_type.isin(@bioactivity_type)")
        # standardize the smiles & clean possible solvents & salts from the string
        .pipe(flag_missing_canonical_smiles)
        .assign(standard_smiles=lambda x: stdzer(x["canonical_smiles"]))
        .pipe(flag_missing_standard_smiles)
        # .dropna(subset=["standard_smiles"])  # drop if no structure is found
        .pipe(flag_salt_or_solvent_removal)
        .assign(final_smiles=lambda x: x["standard_smiles"].apply(clean_mixtures))
        .drop(columns="standard_smiles")
        .rename(columns={"final_smiles": "standard_smiles"})
        .drop(columns=[c for c in cols_to_remove_post_standardization if c in full_df.columns])
        .reset_index(drop=True)
        .copy()
    )
    # TODO: The curate bioactivity errors steps need to be performed here, after the standardization

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
        df = df.assign(
            undefined_stereocenters=lambda x: x["standard_smiles"]
            .apply(find_undefined_stereocenters)
            .apply(len)
        ).pipe(flag_undefined_stereochemistry)
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
    df = flag_inter_document_duplication(df)

    # This part needs to be removed prior to data aggregation. Here, we have either
    # inorganic compounds (SMILES removed from the salt removal step), mixtures (SMILES with "."), or
    # compounds with missing pchembl values, which are needed for the statistics
    missing_smiles_patt = r"Missing Standard SMILES|Missing SMILES|Mixture in SMILES"
    removed_subset = df.query(  # Need SMILES & pchembl_value for aggregation; remove rows with missing them
        "data_dropping_comment.str.contains(@missing_smiles_patt, na=False, regex=True) | "
        "pchembl_value.isna() | "
        r"standard_smiles.str.contains('^\.+$', regex=True)"
    ).copy()
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
    max_assay_match: bool = False,  # This will be driven by perform_assay_match
    output_path: Optional[Union[str, Path]] = None,
    compound_equality: Literal["mixed_fp", "connectivity"] = "connectivity",
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
        max_assay_match: If True, includes assay metadata fields used by Landrum & Riniker, 2024
            for the max assay match. Defaults to False.
        output_path: path to save the aggregated data
        compound_equality: How to identify same compounds in the dataset. If "mixed_fp", uses
            mixed fingerprints (ECFP4 + RDKitFP) to identify same compounds. If "connectivity",
            uses the first part of the InChI key (connectivity) to identify same compounds. Defaults to "connectivity".

    Returns:
        pd.DataFrame: the aggregated data
    """
    current_extra_id_cols = list(extra_id_cols)  # mutable copy

    if max_assay_match:
        logger.info(
            f"Assay metadata matching for aggregation is enabled. Adding fields to ID columns: {DEFAULT_ASSAY_MATCH_FIELDS}"
        )
        # Ensure these columns exist in the DataFrame
        missing_metadata_cols = [col for col in DEFAULT_ASSAY_MATCH_FIELDS if col not in df.columns]
        if missing_metadata_cols:
            logger.warning(
                f"Assay metadata matching enabled for aggregation, but the following "
                f"DEFAULT_ASSAY_MATCH_FIELDS are missing from the DataFrame and will be ignored: {missing_metadata_cols}"
            )
            fields_to_add = [col for col in DEFAULT_ASSAY_MATCH_FIELDS if col in df.columns]
        else:  # Use only the fields that are actually present
            fields_to_add = DEFAULT_ASSAY_MATCH_FIELDS

        current_extra_id_cols.extend(fields_to_add)
        current_extra_id_cols = sorted(list(set(current_extra_id_cols)))
        logger.info(f"ID columns for aggregation: {current_extra_id_cols}")

    connectivity_writer = InchiHandling(convert_to="connectivity", n_jobs=4, progress=True, from_smi=True)

    if compound_equality == "mixed_fp":
        fps = calculate_mixed_FPs(  # Fingerprints are calculated to identify same molecules in the dataset
            df["standard_smiles"].tolist(), n_jobs=4, morgan_kwargs={"useChirality": chirality}
        )
        df = df.assign(id_array=fps)
    elif compound_equality == "connectivity":
        df = df.assign(id_array=lambda x: connectivity_writer(x["standard_smiles"].tolist()))
    # TODO: if other sensible cpd equality choices come along in the future, we can add here...
    else:
        raise ValueError(
            f"Invalid compound_pairing value: {compound_equality}. " "Expected 'mixed_fp' or 'connectivity'."
        )

    # Here we have a repeat index for compounds across all fetched data. Processing which repeats
    # get aggregated (e.g.: same target ID, same `extra_id_cols`, etc) is done in `process_repeat_mols`.
    repeats_idxs = repeated_indices_from_array_series(df["id_array"])

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
        extra_id_cols=current_extra_id_cols,  # Use the potentially extended list
        chirality=chirality,
        extra_multival_cols=include_metadata,
        aggregate_mutants=aggregate_mutants,
    )

    # TODO: we calculate connectivities twice if compound_equality == 'connectivity', which
    # is not ideal. I tried creating a SMILES mapping, but we also canonicalize them inside
    # `process_repeat_mols`, so it doesn't work. Would be nice to fix this in the future
    final_data = final_data.assign(connectivity=lambda x: connectivity_writer(x["smiles"].tolist()))

    # reorder the columns so that connectivity comes first and processing & droppiong comes last
    xtra_cols = [DATA_PROCESSING_COMMENT, DATA_DROPPING_COMMENT]
    cols = ["connectivity", *final_data.columns.difference(["connectivity"] + xtra_cols).tolist(), *xtra_cols]
    final_data = final_data[cols]

    if aggregate_mutants:
        col_subset_dupli_warning = ["connectivity", "target_chembl_id", *extra_id_cols]
    else:
        col_subset_dupli_warning = ["connectivity", "mutation", "target_chembl_id", *extra_id_cols]

    _limit = 15  # limit in the string length for the warning/info logging
    logging_subset = [
        *col_subset_dupli_warning,
        "molecule_chembl_id",
        "assay_chembl_id",
        "pchembl_value_mean",
    ]

    # Based on the ID columns, we shouldn't have any duplicates. This warning is a safeguard
    duplics_for_warning = final_data.duplicated(subset=col_subset_dupli_warning)
    if duplics_for_warning.any():
        dupli_subset = (
            final_data[duplics_for_warning]
            .loc[:, logging_subset]
            .sort_values(
                by=["target_chembl_id", "connectivity", "pchembl_value_mean"],
                ascending=[True, True, False],
            )
        )
        truncated_df = dupli_subset.applymap(
            lambda x: str(x)[:_limit] + "..." if len(str(x)) > _limit else str(x)
        ).head(10)
        logger.warning(
            f"There two or more compounds matching the ID columns {col_subset_dupli_warning} "
            "This is not intentional, please further inspect the collected dataset. Here's a sample "
            "of the repeated entries:\n"
            f"{truncated_df.to_string(index=False)}"
        )

    # Additional safeguard to ensure proper handling of the output by the user prior to modeling
    target_cpd_col_subset = ["target_chembl_id", "connectivity"]
    duplics_for_info = final_data.duplicated(subset=target_cpd_col_subset)
    if duplics_for_info.any():
        dupli_subset = (
            final_data[duplics_for_info]
            .loc[:, logging_subset]
            .sort_values(by=target_cpd_col_subset + ["pchembl_value_mean"], ascending=[True, True, False])
        )
        truncated_df = dupli_subset.applymap(
            lambda x: str(x)[:_limit] + "..." if len(str(x)) > _limit else str(x)
        ).head(10)
        logger.info(
            "There are two or more repeated compound-target readouts (based on `connectivity` & `target_chembl_id`) "
            "without considering other ID columns. This is a result of your aggregation criteria. Make "
            "sure to differ these data points in your modeling pipeline by including information of your other id_columns, "
            "or resolve these compound-target repeats prior to modeling. Here's a sample of the repeated entries:\n"
            f"{truncated_df.to_string(index=False)}"
        )

    if output_path is not None:
        save_dataframe(final_data, output_path)

    return final_data
