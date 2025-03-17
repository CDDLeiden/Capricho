"""Module holding functionalities for the ChEMBL API."""

import re
from typing import List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ..logger import logger
from .api.downloader import get_full_activity_data_sql
from .api.webresource import get_full_activity_data
from .exceptions import BioactivitiesNotFoundError


def convert_to_log10(df):
    """Function to be applied to the whole DataFrame. Will convert the standard_value
    column to pchembl_value column, if the standard_units are in nM, µM or uM.
    Args:
        df: a bioactivity DataFrame. e.g.: output from `get_activity_table`.
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
    column to pchembl_value column if the standard_units are in nM, µM or uM.
    Args:
        bioactivities_df: bioactivity dataframe, e.g.: output from `get_activity_table`.
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


def get_full_activity_workflow(
    molecule_chembl_ids: Optional[Union[list, str]] = None,
    target_chembl_ids: Optional[Union[list, str]] = None,
    assay_chembl_ids: Optional[Union[list, str]] = None,
    document_chembl_ids: Optional[Union[list, str]] = None,
    standard_relation: Optional[List[str]] = None,
    standard_type: Optional[List[str]] = None,
    confidence_scores: Union[list, Tuple] = (9, 8),
    assay_types: Union[list, Tuple] = ("B", "F"),
    chembl_version: Optional[int] = None,
    include_null_values: bool = False,
    additional_fields: Optional[List[str]] = None,
    prefix: Optional[Sequence[str]] = None,
    version: Optional[Union[int, str]] = None,
    calculate_pchembl: bool = False,
    backend: Literal["downloader", "webresource"] = "downloader",
):
    if backend == "downloader":
        bioactivities_df = get_full_activity_data_sql()
    elif backend == "webresource":
        bioactivities_df = get_full_activity_data()

    if bioactivities_df.empty:
        raise BioactivitiesNotFoundError("No bioactivities found for the given query.")
