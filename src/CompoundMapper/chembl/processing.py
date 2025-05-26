"""Module holding functionalities for the ChEMBL API."""

import re
from itertools import combinations
from math import isclose
from typing import List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from ..core.pandas_helper import add_comment
from ..logger import logger
from .api.downloader import get_full_activity_data_sql
from .api.webresource import get_full_activity_data
from .data_flag_functions import (
    flag_calculated_pchembl,
    flag_potential_duplicate,
    flag_with_data_validity_comment,
)
from .exceptions import BioactivitiesNotFoundError


def convert_to_log10(df: pd.DataFrame, save_dropped: bool = False) -> pd.DataFrame:
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
    tmp_df = (
        df.query("standard_units.isin(@desired_units) | standard_units.isna()")
        .copy()
        .pipe(flag_calculated_pchembl)
    )
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
            _info = pchembl_inf_or_nan.loc[pchembl_inf_or_nan.index[:6], debug_cols]
            comment = "Infinite or NaN pchembl_value after calculation"
            if save_dropped:
                logger.info(f"Flagging {len(pchembl_inf_or_nan)} rows: {comment}:\n{_info}")
                tmp_df = add_comment(
                    df=tmp_df,
                    comment=comment,
                    criteria_func=lambda x: x.index.isin(pchembl_inf_or_nan.index),
                    target_column="pchembl_value",  # Dummy target, criteria_func handles selection
                    comment_type="d",
                )
            else:
                logger.info(f"Dropping {len(pchembl_inf_or_nan)} rows: {comment}:\n{_info}")
                tmp_df = tmp_df.drop(pchembl_inf_or_nan.index)
    return tmp_df


# TODO: Consider whether we should do this across pairs identified from fingerprints or molecule IDs
def curate_activity_pairs(
    df: pd.DataFrame,
    mol_id_col: str = "molecule_chembl_id",
    assay_id_col: str = "assay_chembl_id",
    activity_value_col: str = "pchembl_value",
    save_dropped: bool = False,
) -> pd.DataFrame:
    """Curate activity pairs for the same molecule across different assays. Removes pairs
    of measurements if their pChEMBL values differ by 3.0.

    Filter inspired on Landrum & Riniker, 2024, where the authors state:

    > Given the very low probability of two separate experiments producing exactly the same
    results, the exact matches are most likely cases where values from a previous paper are
    copied into a new one; this was discussed in the earlier work by Kramer et al. (10) and
    spot-checked with a number of assay pairs here.

    Args:
        df: DataFrame with bioactivity data.
        mol_id_col: column name for molecule IDs. Defaults to "molecule_chembl_id" for ChEMBL data.
        assay_id_col: column name for assay IDs. Defaults to "assay_chembl_id" for ChEMBL data.
        activity_id_col: column name for activity values. Defaults to "pchembl_value" for ChEMBL data.

    Returns:
        pd.DataFrame: The curated DataFrame.
    """
    if not {mol_id_col, assay_id_col, activity_value_col}.issubset(df.columns):
        logger.warning(
            "Skipping activity pair curation: Required columns "
            f"({mol_id_col}, assay_id_col, 'pchembl_value') not found."
        )
        return df

    rows_to_drop = set()
    # Group by molecule
    for molecule_id, group in df.groupby(mol_id_col):
        # Get unique assays for this molecule
        unique_assays = group[assay_id_col].unique()
        if len(unique_assays) < 2:
            continue  # Not enough assays to form a pair

        # Iterate over all unique pairs of assays for this molecule
        for assay1_id, assay2_id in combinations(unique_assays, 2):
            measurements1 = group[group[assay_id_col] == assay1_id]
            measurements2 = group[group[assay_id_col] == assay2_id]

            # Iterate over all combinations of measurements between the two assays
            for idx1, row1 in measurements1.iterrows():
                for idx2, row2 in measurements2.iterrows():
                    val1 = row1[activity_value_col]
                    val2 = row2[activity_value_col]

                    if pd.isna(val1) or pd.isna(val2):
                        continue

                    # Add a small epsilon for robust floating point comparison >= 3.0
                    if isclose(abs(val1 - val2), 3.0, rel_tol=1e-9, abs_tol=1e-9):
                        rows_to_drop.add(idx1)
                        rows_to_drop.add(idx2)
                        logger.debug(
                            f"Marking rows for molecule {molecule_id}, assays {assay1_id} (value: {val1}) "
                            f"and {assay2_id} (value: {val2}) for removal due to pChEMBL value similarity/difference."
                        )

    if rows_to_drop:
        comment = "potential unit annotation error - pChEMBL values differ by 3.0 for same molecule"
        if save_dropped:
            logger.info(f"Activity Curation: Flagging {len(rows_to_drop)} measurements due to {comment}.")
            # Use a dummy target_column as criteria_func handles selection by index
            df = add_comment(
                df=df,
                comment=comment,
                criteria_func=lambda x: x.index.isin(rows_to_drop),
                target_column=activity_value_col,  # Actual column doesn't matter here
                comment_type="d",
            )
        else:
            logger.info(f"Activity Curation: Removing {len(rows_to_drop)} measurements due to {comment}.")
            df = df.drop(index=list(rows_to_drop))
    return df


def process_bioactivities(
    bioactivities_df: pd.DataFrame,
    calculate_pchembl: bool = True,
    curate_annotation_errors: bool = True,
    require_document_date: bool = False,
    save_dropped: bool = False,
) -> pd.DataFrame:
    """Processes the bioactivities DataFrame. Will convert the standard_value
    column to pchembl_value column if the standard_units are in mM, µM, uM, or nM. If the
    standard_units are in log, the original value in the pchembl_value is preserved.

    Args:
        bioactivities_df: bioactivity dataframe, e.g.: output from `get_activity_table`.
        calculate_pchembl: Whether to calculate pChEMBL values.
        curate_annotation_errors: Whether to apply activity curation based on pChEMBL values
            diverging in exactly 3.0 (indicate possible annotation errors). Defaults to True.
        require_document_date: Whether to filter out activities without a document year.

    Returns:
        pd.DataFrame: the processed bioactivities DataFrame.
    """
    bioactivities_df = (
        bioactivities_df.astype({"standard_value": "float32", "pchembl_value": "float32"})
        .replace({None: np.nan})
        .pipe(flag_with_data_validity_comment)
        # .query("data_validity_comment.isna()")
        .pipe(flag_potential_duplicate)
        # .query("potential_duplicate == 0")
        .drop(columns=["data_validity_comment", "potential_duplicate"])
    )
    with_pchembl = bioactivities_df.query("~pchembl_value.isna()")
    if calculate_pchembl:
        without_pchembl = convert_to_log10(  # if pchembl value not present, calculate it
            bioactivities_df.query("pchembl_value.isna()"), save_dropped=save_dropped
        )
        if not without_pchembl.empty:
            bioactivities_df = pd.concat([with_pchembl, without_pchembl], ignore_index=True)
    else:
        bioactivities_df = with_pchembl

    if curate_annotation_errors:
        bioactivities_df = curate_activity_pairs(bioactivities_df)

    if require_document_date:
        if "year" not in bioactivities_df.columns:
            logger.warning(
                "Document date curation enabled, but 'year' column not found. Skipping this curation."
            )
        else:
            original_count = len(bioactivities_df)
            bioactivities_df = bioactivities_df[bioactivities_df["year"].notna()]
            removed_count = original_count - len(bioactivities_df)
            if removed_count > 0:
                logger.info(
                    f"Document Date Curation: Removed {removed_count} measurements lacking a document year."
                )

    if not save_dropped:
        bioactivities_df = bioactivities_df.assign(  # drop the columns that are flagged
            data_dropping_comment=lambda x: x.data_dropping_comment.replace("", None)
        ).query("data_dropping_comment.isna()")

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
    chembl_release: Optional[int] = None,
    additional_fields: Optional[List[str]] = None,
    prefix: Optional[Sequence[str]] = None,
    version: Optional[Union[int, str]] = None,
    calculate_pchembl: bool = False,
    curate_annotation_errors: bool = True,
    require_document_date: bool = False,
    save_dropped: bool = False,
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
        chembl_release: Not to confuse for `version`. This is the ChEMBL release number used to
            filter the data. Defaults to None.
        additional_fields: `backend=="downloader"` only! "Optional list of additional fields to
            include in the sql query. E.g.: ["vs.sequence"], to retrieve the sequence of the
            variant, if available. Defaults to None.
        prefix: `backend=="downloader"` only! prefix to be used by pystow for storing the data
            on a custom directory. Defaults to None.
        version: `backend=="downloader"` only! version of the ChEMBL database to be downloaded by
            chembl_downloader. If left as None, the latest version will be downloaded. Defaults to None.
        curate_annotation_errors: Whether to apply activity curation based on pChEMBL values diverging
            in exactly 3.0 (indicate possible annotation errors). Defaults to True.
        calculate_pchembl: calculate pChEMBL values for bioactivities reported in nM, µM or uM `standard_unit`
            or -Log|Log `standard_type`. Defaults to False
        save_dropped: If True, rows that would be dropped are kept and flagged with a comment.
            A separate file with these dropped rows might be saved by the calling workflow. Defaults to False.
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
            chembl_release=chembl_release,
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
            chembl_release=chembl_release,
        ).sort_values(by=["molecule_chembl_id", "activity_id", "standard_value"])

    if bioactivities_df.empty:
        raise BioactivitiesNotFoundError(
            "No bioactivities found for the given query: "
            f"molecule_chembl_ids={molecule_chembl_ids}, "
            f"target_chembl_ids={target_chembl_ids}, "
            f"assay_chembl_ids={assay_chembl_ids}, "
            f"document_chembl_ids={document_chembl_ids}, "
        )

    bioactivities_df = process_bioactivities(
        bioactivities_df,
        calculate_pchembl=calculate_pchembl,
        curate_annotation_errors=curate_annotation_errors,
        require_document_date=require_document_date,
        save_dropped=save_dropped,
    )

    if bioactivities_df.empty:
        raise BioactivitiesNotFoundError(
            "No bioactivities found after calculating pchembl_values. To investigate, use "
            "either methods CompoundMapper.chembl.api.downloader.get_full_activity_data_sql or "
            "CompoundMapper.chembl.api.webresource.get_full_activity_data."
        )

    return bioactivities_df
