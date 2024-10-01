"""Module holding functionalities for the ChEMBL API."""

from typing import Optional, Tuple, Union

import numpy as np
import pandas as pd

from ..logger import logger
from .api import (
    assay_info_from_chembl,
    bioactivities_from_chembl,
    get_publications_details,
    molecule_info_from_chembl,
)


def convert_to_log10(df):
    """Function to be applied to the whole DataFrame. Will convert the standard_value
    column to pchembl_value column, if the standard_units are in nM, µM or uM.
    Args:
        df: a bioactivity DataFrame. e.g.: output from `bioactivities_from_chembl`.
    """

    def compute_log(row):
        value = row["standard_value"]
        unit = row["standard_units"]
        if value == 0:  # avoid division by zero
            return np.nan
        if unit == "µM" or unit == "uM":
            value_in_M = value * 1e-6
        elif unit == "nM":
            value_in_M = value * 1e-9
        return -np.log10(value_in_M)

    desired_units = ["nM", "µM", "uM"]
    final_df = []
    for unit in desired_units:  # noqa: B007
        tmp_df = df.query("standard_units == @unit")
        if tmp_df.shape[0] > 0:
            tmp_df = tmp_df.assign(pchembl_value=lambda x: x.apply(compute_log, axis=1))
            pchembl_inf_or_nan = tmp_df.replace([np.inf, -np.inf], np.nan).query("pchembl_value.isna()")
            if not pchembl_inf_or_nan.empty:
                debug_cols = [
                    "target_chembl_id",
                    "assay_chembl_id",
                    "assay_type",
                    "molecule_chembl_id",
                    "standard_units",
                    "standard_value",
                ]
                _info = pchembl_inf_or_nan[debug_cols]
                logger.info(
                    f"Found infinite or NaN values upon calculating pchembl values. Dropping:\n{_info}"
                )
                tmp_df = tmp_df.drop(pchembl_inf_or_nan.index)
            final_df.append(tmp_df)
        else:
            continue
    if final_df == []:
        return None
    else:
        return pd.concat(final_df, ignore_index=True)


def process_bioactivities(bioactivities_df: pd.DataFrame, calculate_pchembl: bool = True) -> pd.DataFrame:
    """Processes the bioactivities DataFrame. Will convert the standard_value
    column to pchembl_value column if the standard_units are in nM, µM or uM.
    Args:
        bioactivities_df: bioactivity dataframe, e.g.: output from `bioactivities_from_chembl`.
    Returns:
        pd.DataFrame: the processed bioactivities DataFrame.
    """
    without_pchembl = None
    bioactivities_df = (
        bioactivities_df.astype({"standard_value": "float32", "pchembl_value": "float32"})
        .replace({None: np.nan})
        .query("data_validity_description.isna()")
        .query("potential_duplicate == 0")
        .drop(columns=["data_validity_description", "potential_duplicate"])
    )
    with_pchembl = bioactivities_df.query("~pchembl_value.isna()")
    if calculate_pchembl:
        without_pchembl = convert_to_log10(  # if pchembl value not present, calculate it
            bioactivities_df.query("pchembl_value.isna()")
        )
        if without_pchembl is not None:
            bioactivities_df = pd.concat([with_pchembl, without_pchembl], ignore_index=True)
    else:
        bioactivities_df = with_pchembl
    return bioactivities_df.reset_index(drop=True)


def fetch_and_filter_workflow(
    molecule_chembl_ids: Optional[list] = None,
    target_chembl_ids: Optional[list] = None,
    assay_chembl_ids: Optional[list] = None,
    document_chembl_ids: Optional[list] = None,
    confidence_scores: Union[list, Tuple] = (9, 8),
    assay_types: Union[list, Tuple] = ("B", "F"),
    calculate_pchembl: bool = True,
    chembl_version: Optional[int] = None,
) -> pd.DataFrame:
    """
    This function retrieves and merges data from ChEMBL for given molecule or target IDs.
    The steps taken to filter this data are:

    1. Fetch bioactivities for the given target or molecule IDs, considering designated
        confidence scores and assay types (binding or functional) using `new_client.activity`
        from them `chembl_webresource_client` package.
    2. If the `calculate_pchembl` parameter is true, use obtained bioactivity points with values
        reported in (standard units) nM, µM or uM to calculate the pChEMBL value.
    3. Extract unique assay IDs from the bioactivities DataFrame & add this information to the
        final DataFrame.

    Args:
        molecule_chembl_ids: list of ChEMBL molecule IDs to fetch data for. Defaults to None.
        target_chembl_ids: list of ChEMBL target IDs to fetch data for. Defaults to None.
        assay_chembl_ids: list of ChEMBL assay IDs to fetch data for. Defaults to None.
        document_chembl_ids: list of ChEMBL document IDs to fetch data for. Defaults to None.
        confidence_scores: list of confidence scores to filter the fetched assay data.
            Defaults to (9, 8).
        assay_types: list of assay types to be fetched from ChEMBL. Defaults to binding (B) and
            functional (F) data.
        calculate_pchembl: calculate pChEMBL values for the bioactivities when it's not present
            for bioactivities reported in nM, µM or uM. Defaults to True.
        chembl_version: specify latest ChEMBL release to extract data from (e.g., 28). Defaults to None.
    Returns:
        pd.DataFrame: Merged DataFrame with molecule, bioactivity, and assay information.
    """
    bioactivities_df = process_bioactivities(
        bioactivities_from_chembl(
            molecule_chembl_ids=molecule_chembl_ids,
            target_chembl_ids=target_chembl_ids,
            assay_chembl_ids=assay_chembl_ids,
            document_chembl_ids=document_chembl_ids,
            standard_relation="=",
            assay_type__in=assay_types,
        ),
        calculate_pchembl=calculate_pchembl,
    )
    unique_assay_ids = bioactivities_df["assay_chembl_id"].unique().tolist()

    logger.debug(f"Fetching assays of type {assay_types} with confidence scores: {list(confidence_scores)}")
    assays_df = assay_info_from_chembl(
        unique_assay_ids, confidence_scores=list(confidence_scores), assay_type__in=assay_types
    )

    merged_df = pd.merge(
        bioactivities_df,
        assays_df,
        on=["assay_chembl_id", "target_chembl_id", "assay_type"],
        how="right",  # keep only the bioactivities with respective assays
    ).drop_duplicates()
    logger.debug(f"Columns in the merged DataFrame: {merged_df.columns}")

    # Get the molecule structures & merge to the dataset
    unique_mol_chembl_ids = merged_df["molecule_chembl_id"].unique().tolist()
    logger.debug(f"Fetching structural information for {len(unique_mol_chembl_ids)} molecule ChEMBL IDs.")
    mol_data = molecule_info_from_chembl(unique_mol_chembl_ids)

    # Merge using molecule_chembl_id
    full_df = merged_df.merge(mol_data, on="molecule_chembl_id", how="inner")
    full_df = full_df[  # reorder columns with molecule_chembl_id and canonical_smiles first
        ["molecule_chembl_id", "canonical_smiles"]
        + [col for col in full_df.columns if col not in ["molecule_chembl_id", "canonical_smiles"]]
    ]
    logger.debug(f"Columns in the DataFrame with molecular structures: {full_df.columns}")
    if chembl_version is not None:
        document_ids = full_df["document_chembl_id"].unique().tolist()
        logger.info("Fetching publication details for the documents.")
        publications_df = get_publications_details(document_ids, chembl_version)
        full_df = pd.merge(full_df, publications_df, on="document_chembl_id", how="right")

    return full_df
