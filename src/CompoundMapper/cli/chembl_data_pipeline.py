from inspect import signature
from pathlib import Path
from typing import Literal, Optional, Union

import pandas as pd
from chemFilters.chem.standardizers import ChemStandardizer

from ..chembl.data_flag_functions import (
    flag_duplication,
    flag_insufficient_assay_overlap,
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
    DEFAULT_ASSAY_MATCH_FIELDS,
    MOLECULE_ID,
    TARGET_ID,
)
from ..core.fp_utils import calculate_mixed_FPs
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import process_repeat_mols, repeated_indices_from_array_series
from ..core.stereo import find_undefined_stereocenters
from ..logger import logger


def get_standardize_and_clean_workflow(
    molecule_ids: list[str],
    target_ids: list[str],
    assay_ids: list[str],
    document_ids: list[str],
    calculate_pchembl: bool,
    output_path: Optional[Union[str, Path]],
    confidence_scores: list[str],
    bioactivity_type: list[str],
    standard_relation: list[str],  # TODO: later support data with  >, <, >=, <=...
    assay_types: list[str],
    chembl_release: Optional[int],
    save_not_aggregated: bool = True,
    save_dropped: bool = False,
    save_duplicated: bool = False,
    drop_unassigned_chiral: bool = False,
    version: Optional[Union[int, str]] = None,
    backend: Literal["downloader", "webresource"] = "downloader",
    curate_annotation_errors: bool = True,
    require_doc_date: bool = False,
    max_assay_size: Optional[int] = None,
    min_assay_overlap: int = 0,
    strict_mutant_removal: bool = False,
    keep_flagged_data: bool = False,
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
        output_path: path to save the resulting csv file
        confidence_scores: list of confidence scores (assay-related) to filter data from
        bioactivity_type: list of bioactivity types (assay-related) to filter data from
        standard_relation: standard relation to filter data from. Currently only supports "="
            Defaults to "=".
        chembl_release: latest ChEMBL release to retrieve data from
        save_not_aggregated: whether to save the resulting data to the csv (output_path) before
        save_dropped: whether to save a separate dataframe containing rows that were flagged for dropping.
        save_duplicated: whether to save the duplicated data (if any) to a separate csv file
        drop_unassigned_chiral: whether to drop data points with undefined stereocenters. Defaults to False.
        version: `backend=="downloader"` only! version of the ChEMBL database to be downloaded by
            chembl_downloader. If left as None, the latest version will be downloaded. Defaults to None.
        backend: the backend to be used for fetching the data. If downloader, the ChEMBL sql database
            is downloaded and extracted first. Defaults to "downloader".
        curate_annotation_errors: Whether to apply activity curation based on pChEMBL values diverging
            in exactly 3.0 (indicate possible annotation errors). Defaults to True.
        require_doc_date: Whether to filter out activities without a document year.
        max_assay_size: Maximum number of compounds in an assay. Assays exceeding this size will
            have their activities flagged for removal. Defaults to None (no filtering).
        min_assay_overlap: Minimum number of overlapping compounds between two assays for the same target
                for their activities to be considered. Defaults to None (no filtering).
            strict_mutant_removal: If True, assays with 'mutant', 'mutation', or 'variant' in their
                description will be flagged for removal. Defaults to False.
        keep_flagged_data: If True, data points flagged for dropping (due to various quality
            checks) will be retained in the main DataFrame. The `data_dropping_comment` column
            will still be populated. A warning will be logged. Defaults to False.

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
        save_dropped=save_dropped,
        version=version,
        backend=backend,
    )
    # Filter assays by size
    if max_assay_size is not None and not full_df.empty:
        logger.info(f"Filtering assays with more than {max_assay_size} compounds.")
        # Ensure the necessary columns exist before trying to group or filter
        if ASSAY_ID in full_df.columns and MOLECULE_ID in full_df.columns:
            assay_counts = full_df.groupby(ASSAY_ID)[MOLECULE_ID].nunique()
            assays_to_filter = assay_counts[assay_counts > max_assay_size].index.tolist()

            if assays_to_filter:
                filter_mask = full_df[ASSAY_ID].isin(assays_to_filter)
                num_activities_flagged = filter_mask.sum()
                num_assays_flagged = len(assays_to_filter)

                logger.info(
                    f"Flagging {num_activities_flagged} activities from {num_assays_flagged} assays "
                    f"exceeding max size of {max_assay_size}."
                )
                # Ensure DATA_DROPPING_COMMENT column exists
                if DATA_DROPPING_COMMENT not in full_df.columns:
                    full_df[DATA_DROPPING_COMMENT] = pd.Series(dtype="object")

                comment = f"Assay size exceeds maximum (N > {max_assay_size})"
                existing_comments = full_df.loc[filter_mask, DATA_DROPPING_COMMENT].fillna("")
                new_comments = existing_comments.apply(lambda x: f"{x} & {comment}" if x else comment)
                full_df.loc[filter_mask, DATA_DROPPING_COMMENT] = new_comments
            else:
                logger.info("No assays found exceeding the maximum size.")
        else:
            logger.warning(
                f"Could not filter by assay size. Required columns "
                f"'{ASSAY_ID}' or '{MOLECULE_ID}' not found in DataFrame."
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
            logger=logger,
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

    stdzer = ChemStandardizer(from_smi=True, n_jobs=8, verbose=False)
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

    # Raise error if df is empty after critical processing steps AND we are not saving dropped rows
    if not save_dropped and df.empty:
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
    if drop_unassigned_chiral:
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

    # Documents - since much of ChEMBL is take from medchem literature, if there are two
    # assays from the same document that have the same readout and are measured against
    # the same target, we can be pretty sure tha they are *not* equivalent;
    df_size = df.shape[0]
    col_subset = [  # Drop different assays that have the same exact pchembl value - probably duplicate
        "molecule_chembl_id",
        "standard_smiles",
        "canonical_smiles",
        "pchembl_value",
        "standard_relation",
        "target_chembl_id",
        "target_organism",
    ]
    # Handle duplicated data
    duplicated = df.duplicated(subset=col_subset, keep=False)
    df = flag_duplication(df, dupli_id_subset=col_subset)
    if duplicated.any():
        df.drop_duplicates(subset=col_subset, keep="first", inplace=True)
        if save_duplicated and output_path is not None:
            df[duplicated].to_csv(output_path.with_stem(f"{output_path.stem}_duplicated"), index=False)
        if df_size != df.shape[0]:
            logger.info(f"Dropped {df_size - df.shape[0]} duplicates.")

    if save_dropped and output_path is not None:
        df.assign(data_dropping_comment=lambda x: x.data_dropping_comment.replace("", None)).query(
            "~data_dropping_comment.isna()"
        ).to_csv(output_path.with_stem(f"{output_path.stem}_removed"), index=False)

    if keep_flagged_data:
        logger.warning(
            "Retaining flagged data points in the main dataset. "
            "This dataset includes entries that would normally be dropped due to quality flags "
            "and is generally not recommended. The 'data_dropping_comment' column indicates the reasons."
        )
        missing_smiles_patt = "Missing (standard )SMILES|Missing SMILES"
        df = df.query(  # Need SMILES & pchembl_value for aggregation; remove rows with missing them
            "~data_dropping_comment.str.contains(@missing_smiles_patt, na=False, regex=True)"
        ).query("~pchembl_value.isna()")

    else:
        df = df.assign(data_dropping_comment=lambda x: x.data_dropping_comment.replace("", None)).query(
            "data_dropping_comment.isna()"
        )

    if save_not_aggregated and output_path is not None:
        df.to_csv(output_path.with_stem(f"{output_path.stem}_not_aggregated"), index=False)

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

    Returns:
        pd.DataFrame: the aggregated data
    """
    current_extra_id_cols = list(extra_id_cols)  # Make a mutable copy

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
        logger.info(f"Effective ID columns for aggregation: {current_extra_id_cols}")

    fps = calculate_mixed_FPs(  # Fingerprints are calculated to identify same molecules in the dataset
        df["standard_smiles"].tolist(), n_jobs=4, morgan_kwargs={"useChirality": chirality}
    )
    df = df.assign(fps=fps)
    repeats_idxs = repeated_indices_from_array_series(df["fps"])

    include_metadata = ["doc_type", "doi", "journal", "year", "chembl_release", *extra_multival_cols]

    final_data = process_repeat_mols(
        df,
        repeats_idxs,
        solve_strat="keep",
        extra_id_cols=current_extra_id_cols,  # Use the potentially extended list
        chirality=chirality,
        extra_multival_cols=include_metadata,
        aggregate_mutants=aggregate_mutants,
    )

    if output_path is not None:
        final_data.to_csv(output_path, index=False)

    return final_data
