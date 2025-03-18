"""Module holding functionalities for the ChEMBL API."""

import re
from typing import List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ..logger import logger
from .api.downloader import get_full_activity_data_sql
from .api.webresource import get_full_activity_data
from .exceptions import BioactivitiesNotFoundError


def convert_to_log10(df) -> pd.DataFrame:
    """Function to be applied to the whole DataFrame. Will convert the standard_value
    column to pchembl_value column, if the standard_units are in nM, µM or uM.

    Args:
        df: a bioactivity DataFrame. e.g.: output from `get_activity_table`.

    Returns:
        pd.DataFrame: the DataFrame with the pchembl_value column added. If the data didn't
            meet the pchembl_value calculation criteria (standard_units not in nM, µM or uM) |
            (standard_type.str.contains(log)), it can return an empty DataFrame.
    """

    def compute_log(row):
        # Sometimes standard_type is already log transformed;
        if re.findall("log", row["standard_type"], re.IGNORECASE):
            if row["standard_value"] < 0:
                return -row["standard_value"]  # log10 transformed, but not negative yet ?
            else:
                return row["standard_value"]

        unit = row["standard_units"]
        if pd.isna(unit):
            return np.nan

        value = row["standard_value"]
        if value == 0:  # avoid division by zero
            return np.nan
        if unit == "mM":
            value_in_M = value * 1e-3
        elif unit == "µM" or unit == "uM":
            value_in_M = value * 1e-6
        elif unit == "nM":
            value_in_M = value * 1e-9

        return -np.log10(value_in_M)

    desired_units = ["nM", "µM", "uM", "mM"]  # noqa: F841
    # allow na values in case value is log transformed
    tmp_df = df.query("standard_units.isin(@desired_units) | standard_units.isna()").copy()
    if "pchembl_value" not in tmp_df.columns:
        raise ValueError("pchembl_value column not found in input DataFrame.")
    elif (~tmp_df.pchembl_value.isna()).any():
        raise ValueError(
            "input DataFrame should only have pchembl_value.isna() values. If not, use "
            "chembl.processing.process_bioactivities instead."
        )
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
            logger.info(f"Found infinite or NaN values upon calculating pchembl values. Dropping:\n{_info}")
            tmp_df = tmp_df.drop(pchembl_inf_or_nan.index)
    return tmp_df


def process_bioactivities(bioactivities_df: pd.DataFrame, calculate_pchembl: bool = True) -> pd.DataFrame:
    """Processes the bioactivities DataFrame. Will convert the standard_value
    column to pchembl_value column if the standard_units are in mM, µM, uM, or nM. If the
    standard_units are in log, the original value in the pchembl_value is preserved.

    Args:
        bioactivities_df: bioactivity dataframe, e.g.: output from `get_activity_table`.

    Returns:
        pd.DataFrame: the processed bioactivities DataFrame.
    """
    print(bioactivities_df.columns)
    bioactivities_df = (
        bioactivities_df.astype({"standard_value": "float32", "pchembl_value": "float32"})
        .replace({None: np.nan})
        .query("data_validity_comment.isna()")
        .query("potential_duplicate == 0")
        .drop(columns=["data_validity_comment", "potential_duplicate"])
    )
    with_pchembl = bioactivities_df.query("~pchembl_value.isna()")
    if calculate_pchembl:
        without_pchembl = convert_to_log10(  # if pchembl value not present, calculate it
            bioactivities_df.query("pchembl_value.isna()")
        )
        if not without_pchembl.empty:
            bioactivities_df = pd.concat([with_pchembl, without_pchembl], ignore_index=True)
    else:
        bioactivities_df = with_pchembl
    return bioactivities_df.reset_index(drop=True)


def get_bioactivities_workflow(
    molecule_chembl_ids: Optional[Union[list, str]] = None,
    target_chembl_ids: Optional[Union[list, str]] = None,
    assay_chembl_ids: Optional[Union[list, str]] = None,
    document_chembl_ids: Optional[Union[list, str]] = None,
    standard_relation: Optional[List[str]] = None,
    standard_type: Optional[List[str]] = None,
    confidence_scores: Union[list, Tuple] = (9, 8),
    assay_types: Union[list, Tuple] = ("B", "F"),
    chembl_version: Optional[int] = None,
    additional_fields: Optional[List[str]] = None,
    prefix: Optional[Sequence[str]] = None,
    version: Optional[Union[int, str]] = None,
    calculate_pchembl: bool = False,
    backend: Literal["downloader", "webresource"] = "downloader",
):
    """Perform the first step of the bioactivity data retrieval workflow. These are:

    1. Get ChEMBL data using any of the input identifiers: molecule_chembl_ids, target_chembl_ids,
    assay_chembl_ids, or document_chembl_ids. Additional filters are supported by other input
    parameters.

    2. Once the data is retrieved, the bioactivities are processed according to the `calculate_pchembl`
    parameter. ChEMBL calculates pChEMBL values for activity data with the following criteria:

    - "standard_type" in "[IC50", "XC50", "EC50", "AC50", "Ki", "Kd", "Potency", "ED50"];
    - "standard_relation" == "=" & "standard_units" == "nM";
    - "standard_value" > 0 & "data_validity_comment"].isnull()) | "data_validity_comment" == "Manually validated"

    By passing this parameter to True, pchembl values will be calculated for bioactivities reported in
    nM, µM or uM `standard_unit`, or -Log|Log `standard_type`.

    3. The first quality filter is applied. Data containing "data_validity_description" or
    "potential_duplicate" flags are immediatelly removed from the DataFrame.

    Args:
        molecule_chembl_ids: list of ChEMBL molecule IDs to fetch data for. Defaults to None.
        target_chembl_ids: list of ChEMBL target IDs to fetch data for. Defaults to None.
        assay_chembl_ids: list of ChEMBL assay IDs to fetch data for. Defaults to None.
        document_chembl_ids: list of ChEMBL document IDs to fetch data for. Defaults to None.
        standard_relation: Optional filter for standard relation types (e.g., ["=", "<", ">"])
        standard_type: Optional filter for activity types (e.g., ["IC50", "Ki", "EC50"])
        confidence_scores: list of confidence scores to filter the fetched assay data.
            Defaults to (9, 8).
        assay_types: list of assay types to be fetched from ChEMBL. Defaults to binding (B) and
            functional (F) data.
        chembl_version: Not to confuse for `version`. This is the ChEMBL release number used to
            filter the data. Defaults to None.
        additional_fields: `backend=="downloader"` only! "Optional list of additional fields to
            include in the sql query. E.g.: ["vs.sequence"], to retrieve the sequence of the
            variant, if available. Defaults to None.
        prefix: `backend=="downloader"` only! prefix to be used by pystow for storing the data
            on a custom directory. Defaults to None.
        version: `backend=="downloader"` only! version of the ChEMBL database to be downloaded by
            chembl_downloader. If left as None, the latest version will be downloaded. Defaults to None.
        calculate_pchembl: calculate pChEMBL values for bioactivities reported in nM, µM or uM `standard_unit`
            or -Log|Log `standard_type`. Defaults to False
        backend: the backend to be used for fetching the data. If downloader, the ChEMBL sql database
            is downloaded and extracted first. Defaults to "downloader".

    Raises:
        BioactivitiesNotFoundError: If the retrieved bioactivity dataframe is empty.
    """
    if backend == "downloader":
        bioactivities_df = get_full_activity_data_sql(
            molecule_chembl_ids=molecule_chembl_ids,
            target_chembl_ids=target_chembl_ids,
            assay_chembl_ids=assay_chembl_ids,
            document_chembl_ids=document_chembl_ids,
            standard_relation=standard_relation,
            standard_type=standard_type,
            confidence_scores=confidence_scores,
            assay_types=assay_types,
            chembl_version=chembl_version,
            additional_fields=additional_fields,
            prefix=prefix,
            version=version,
        )
    elif backend == "webresource":
        bioactivities_df = get_full_activity_data(
            molecule_chembl_ids=molecule_chembl_ids,
            target_chembl_ids=target_chembl_ids,
            assay_chembl_ids=assay_chembl_ids,
            document_chembl_ids=document_chembl_ids,
            confidence_scores=confidence_scores,
            assay_types=assay_types,
            chembl_version=chembl_version,
        ).sort_values(by=["molecule_chembl_id", "activity_id", "standard_value"])

    if bioactivities_df.empty:
        raise BioactivitiesNotFoundError(
            "No bioactivities found for the given query: "
            f"molecule_chembl_ids={molecule_chembl_ids}, "
            f"target_chembl_ids={target_chembl_ids}, "
            f"assay_chembl_ids={assay_chembl_ids}, "
            f"document_chembl_ids={document_chembl_ids}, "
        )

    bioactivities_df = process_bioactivities(bioactivities_df, calculate_pchembl=calculate_pchembl)

    if bioactivities_df.empty:
        raise BioactivitiesNotFoundError(
            "No bioactivities found after calculating pchembl_values. To investigate, use "
            "either methods CompoundMapper.chembl.api.downloader.get_full_activity_data_sql or "
            "CompoundMapper.chembl.api.webresource.get_full_activity_data."
        )

    return bioactivities_df
