from inspect import signature
from pathlib import Path
from typing import Optional, Union

import numpy as np
from chemFilters.chem.standardizers import ChemStandardizer

from ..chembl.exceptions import BioactivitiesNotFoundError
from ..chembl.processing import fetch_and_filter_workflow
from ..core.fp_utils import calculate_mixed_FPs
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import process_repeat_mols, repeated_indices_from_array_series
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
    chembl_version: int,
    save_not_aggregated: bool,
    save_duplicated: bool = False,
    add_document_info: bool = True,
) -> None:
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
        chembl_version: latest ChEMBL release to retrieve data from
        save_not_aggregated: whether to save the resulting data to the csv (output_path) before
        save_duplicated: whether to save the duplicated data (if any) to a separate csv file
        add_document_info: whether to add publication-related fields to the final DataFrame. Setting
            to True, will require one less query to be made to ChEMBL, but fields like `year` will be
            lacking. Defaults to True.

    Returns:
        pd.DataFrame: the filtered, standardized, and cleaned data
    """
    if output_path is not None:
        if isinstance(output_path, str):
            output_path = Path(output_path)

    # since we work with pchembl values, standard values reported as -pXC50, -Log XC50, etc. will be renamed
    bioactivity_type_rename_dict = {
        **{f"p{bio}": bio for bio in bioactivity_type},
        **{f"Log {bio}": bio for bio in bioactivity_type},
        **{f"-Log {bio}": bio for bio in bioactivity_type},
    }

    full_df = fetch_and_filter_workflow(
        molecule_chembl_ids=molecule_ids,
        target_chembl_ids=target_ids,
        assay_chembl_ids=assay_ids,
        document_chembl_ids=document_ids,
        confidence_scores=confidence_scores,
        assay_types=assay_types,
        calculate_pchembl=calculate_pchembl,
        chembl_version=chembl_version,
        add_document_info=add_document_info,
    )

    # drop rows without chemical structures
    no_smiles_mask = full_df.canonical_smiles.isna()
    if no_smiles_mask.any():
        _info = full_df[no_smiles_mask].iloc[:, :6]
        logger.info(f"Dropping rows with missing canonical smiles:\n{_info}")
        full_df = full_df.drop(index=_info.index).reset_index(drop=True)

    stdzer = ChemStandardizer(from_smi=True, n_jobs=8, verbose=False)
    queried_df = (
        full_df.assign(  # rename bioactivities & filter by preferred bioactivity type
            standard_type=lambda x: x["standard_type"].replace(bioactivity_type_rename_dict)
        )
        .query("standard_type.isin(@bioactivity_type)")
        # standardize the smiles & clean possible solvents & salts from the string
        .assign(standard_smiles=lambda x: stdzer(x["canonical_smiles"]))
        .dropna(subset=["standard_smiles"])  # drop if no structure is found
        .query("standard_smiles.notna()")
        .assign(final_smiles=lambda x: x["standard_smiles"].apply(clean_mixtures))
        .drop(columns="standard_smiles")
        .rename(columns={"final_smiles": "standard_smiles"})
        # Drop cols that won't be used; we'll use `standard_<colname>` instead
        .drop(
            columns=[
                "type",
                "relation",
                "units",
                "value",
                "standard_value",  # we'll use pchembl instead
                "type",
                # "description",
            ]
        )
        .reset_index(drop=True)
        .copy()
    )
    # drop rows with missing pchembl values
    no_pchembl_idxs = queried_df.query("pchembl_value.isna()").index
    logger.info(f"Dropping {len(no_pchembl_idxs)} rows missing pchembl values.")
    queried_df = queried_df.drop(index=no_pchembl_idxs).reset_index(drop=True)

    if queried_df.empty:
        func_params = signature(get_standardize_and_clean_workflow).parameters
        local_vars = locals()
        parameters = {name: local_vars[name] for name in func_params}
        raise BioactivitiesNotFoundError(parameters=parameters)

    # find remaining mixtures in the data
    mask = queried_df["standard_smiles"].str.contains(".", regex=False)
    n_mixtures = mask.sum()
    if n_mixtures > 0:
        logger.info(f"Number of mixtures: {mask.sum()}")
        mixture_idxs = np.where(mask)[0]  # drop where smiles contain mixtures
        queried_df = queried_df.drop(index=mixture_idxs).reset_index(drop=True)

    # Handle duplicated data
    df_size = queried_df.shape[0]
    col_subset = [  # Drop different assays that have the same exact pchembl value - probably duplicate
        "molecule_chembl_id",
        "standard_smiles",
        "canonical_smiles",
        "pchembl_value",
        "standard_relation",
        "target_chembl_id",
        "target_organism",
    ]
    duplicated = queried_df.duplicated(subset=col_subset, keep=False)
    if duplicated.any():
        queried_df.drop_duplicates(
            subset=col_subset,
            keep="first",
            inplace=True,
        )
        if save_duplicated and output_path is not None:
            queried_df[duplicated].to_csv(
                output_path.with_stem(f"{output_path.stem}_duplicated"), index=False
            )
        if df_size != queried_df.shape[0]:
            logger.info(f"Dropped {df_size - queried_df.shape[0]} duplicates.")

    # TODO: in the future, implement the following performed in Greg Landrum's max curation paper:
    # Drop both the assays where we have exactly 3 pchembl value difference - somebody probably did a uM / nM confusion

    # Documents - since much of ChEMBL is take from medchem literature, if there are two
    # assays from the same document that have the same readout and are measured against
    # the same target, we can be pretty sure tha they are *not* equivalent;

    if save_not_aggregated and output_path is not None:
        queried_df.to_csv(output_path.with_stem(f"{output_path.stem}_not_aggregated"), index=False)

    return queried_df


def aggregate_data(
    df,
    chirality: bool,
    chembl_version: int,
    metadata_cols: list[str],
    extra_id_cols: list[str],
    aggregate_mutants: bool = False,
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
        chembl_version: latest ChEMBL release to retrieve data from
        metadata_cols: additional metadata columns to keep in the final dataframe. Metadata will
            be saved separated by a semicolon whenever aggregation is performed.
        extra_id_cols: additional columns to use as identifiers for the aggregation. Passing
            `["assay_chembl_id"]` to this argument, for example, will only aggregate the data
            if the compound is the same and the assay is the same.
        aggregate_mutants: if true, will aggregate data solely based on the target_chembl_id,
            regardless of the variant_sequence flag in ChEMBL. Defaults to False.
        output_path: path to save the aggregated data

    Returns:
        pd.DataFrame: the aggregated data
    """
    fps = calculate_mixed_FPs(  # Fingerprints are calculated to identify same molecules in the dataset
        df["standard_smiles"].tolist(), n_jobs=4, morgan_kwargs={"useChirality": chirality}
    )
    df = df.assign(fps=fps)
    repeats_idxs = repeated_indices_from_array_series(df["fps"])

    if chembl_version is not None:
        include_metadata = ["doc_type", "doi", "journal", "year", "chembl_release", *metadata_cols]
    else:
        include_metadata = metadata_cols

    final_data = process_repeat_mols(
        df,
        repeats_idxs,
        solve_strat="keep",
        extra_id_cols=extra_id_cols,
        chirality=chirality,
        extra_multival_cols=include_metadata,
        aggregate_mutants=aggregate_mutants,
    )

    if output_path is not None:
        final_data.to_csv(output_path, index=False)

    return final_data
