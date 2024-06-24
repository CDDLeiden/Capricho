"""Module containing helper functions for processing repeated elements in a DataFrame"""

from typing import Iterable, List

import numpy as np
import pandas as pd

from ..logger import logger
from .pandas_helper import aggr_val_series, apply_func_grpd, assign_stats


def repeated_indices_from_IDs_df(df: pd.DataFrame, columns: list) -> List[List[int]]:
    """Find repeated indices for given columns in a DataFrame.

    Args:
        df (pd.DataFrame): The DataFrame to search for repeats.
        columns (list): List of column names with external (non-chembl) IDs to identify
            repeats across.; E.g.: ["JUMP_ID", "target_chembl_id"]

    Returns:
        list: A list of lists containing indices of repeated rows based on specified columns.
    """
    # Validate input
    if not all(col in df for col in columns):
        raise ValueError("One or more specified columns do not exist in the DataFrame.")

    # Concatenate values of specified columns into a single series
    concatenated_series = df[columns].astype(str).agg("-".join, axis=1)

    # Find duplicates using numpy operations
    def find_duplicate_index(series: pd.Series) -> List[List[int]]:
        arr = series.to_numpy()
        sidx = np.lexsort(arr.reshape(1, -1))
        sorted_series = series.iloc[sidx]
        sorted_arr = sorted_series.to_numpy()
        # Identify duplicates
        duplicates_mask = np.concatenate(([False], sorted_arr[1:] == sorted_arr[:-1], [False]))
        idx = np.flatnonzero(duplicates_mask[1:] != duplicates_mask[:-1])
        # Extract original indices for duplicates
        sorted_indices = sorted_series.index.tolist()
        return [sorted_indices[i:j] for i, j in zip(idx[::2], idx[1::2] + 1)]

    final_repeat_idxs = find_duplicate_index(concatenated_series)
    return final_repeat_idxs


def repeated_indices_from_array_series(series: Iterable) -> List[List[int]]:
    """Function to find repeated arrays from a list of arrays"""

    def find_duplicate_index(series: np.array) -> List[List[int]]:
        """Group indices of duplicate rows
        From https://stackoverflow.com/a/46629623

        Args:
            arr (numpy array): array of fingerprints

        Returns:
            list[list[int]]: list of lists of indices of duplicate rows
        """
        # Sort by rows
        arr = np.vstack(series)
        sidx = np.lexsort(arr.T)
        b = arr[sidx]
        # Get unique row mask
        m = np.concatenate(([False], (b[1:] == b[:-1]).all(1), [False]))
        # Get start and stop indices for each group of duplicates
        idx = np.flatnonzero(m[1:] != m[:-1])
        # Get sorted indices
        sort_idxs = series.index[sidx].tolist()
        # Return list of lists of indices of duplicate rows
        return [sort_idxs[i:j] for i, j in zip(idx[::2], idx[1::2] + 1)]

    final_repeat_idxs = find_duplicate_index(series)
    return final_repeat_idxs


def process_repeat_mols(
    df: pd.DataFrame,
    repeat_element_idxs: List[List[int]],
    solve_strat: str = "keep",
    extra_id_cols: List[str] = [],
) -> pd.DataFrame:
    """Process the dataframe according to repeated elements identified
    with the function `find_repeated_arr_from_series`. The standard criteria here
    will be that molecules with the same Fingerprint representation will be treated
    as a single entity, and will have their values aggregated. Upon aggregation, if the
    min & max values differ 1 or more log units, then those samples will be remioved
    from the dataset. Otherwise, values will be aggregated and a new column will be
    assigned, called `might_rancemic`. This column will be a boolean, indicating
    whether the molecule might be rancemic or not.

    Args:
        df: dataframe with the bioactivity data
        repeat_element_idxs: list of indices of repeated elements in the dataframe.
        solve_strat: strategy to solve the repeated elements. If 'drop', then both the
            points within >= 1 log unit difference will be dropped. If 'keep', then
            no values will be dropped.
        extra_id_cols: list of extra identification columns you might have for your own
            compounds that you'd like to use to avoid mixing data & to keep in the final
            dataframe.

    Returns:
        df: dataframe with the repeated elements processed.
    """
    df = df.copy()
    repeat_mapping = {}
    for idx in range(len(repeat_element_idxs)):
        for i in repeat_element_idxs[idx]:
            repeat_mapping[i] = idx
    df = df.assign(repeat_mapping=lambda x: x.index.map(repeat_mapping))
    repeat_subset = df.query("~repeat_mapping.isna()").assign(
        pchembl_value=lambda x: x.pchembl_value.astype(str)
    )
    numeric_activity = (
        # concatenate grouped values and convert to numeric arrays
        repeat_subset.groupby(["repeat_mapping"])["pchembl_value"]
        .apply(lambda x: ";".join(x))
        .str.split(";")
        .apply(lambda x: np.array(x).astype(float))
    )
    max_series = numeric_activity.apply(lambda x: np.max(x))
    min_series = numeric_activity.apply(lambda x: np.min(x))
    distance_series = max_series - min_series

    # Will drop the repeats with more than 1 log unit difference
    high_diff_repeats = np.where(distance_series >= 1)[0]
    points_dropped = len(repeat_subset["repeat_mapping"].isin(high_diff_repeats))
    logger.info(f"Found {len(high_diff_repeats)} repeats with more than 1 log unit difference.")
    logger.info(f"Maximum difference between min & max values: {np.max(distance_series)}")
    if solve_strat == "drop":
        logger.info(f"{points_dropped} points will be removed from the dataset")

    id_cols = [*extra_id_cols, "repeat_mapping", "target_chembl_id"]
    multival_cols = [
        "standard_smiles",
        "pchembl_value",
        "assay_chembl_id",
        "assay_description",
        "activity_id",
        "assay_type",
        "standard_type",
        "confidence_score",
        "standard_relation",
        "target_organism",
        "molecule_chembl_id",
        "document_chembl_id",
        "assay_tissue",
        "assay_cell_type",
        "relationship_description",
        "variant_sequence",
        "indication_class",
        "max_phase",
        "oral",
        "prodrug",
        "withdrawn_flag",
    ]
    final_cols = id_cols + multival_cols
    grouped = repeat_subset.groupby(id_cols)
    updated_vals = apply_func_grpd(grouped, aggr_val_series, id_cols, *multival_cols)
    updated_df = assign_stats(updated_vals, value_col="pchembl_value").merge(df, on=id_cols, how="left")
    rename_cols = {c: c.rstrip("_x") for c in updated_df.columns if c.endswith("_x")}
    todrop_cols = [c for c in updated_df.columns if c.endswith("_y")]
    updated_df = updated_df.drop(columns=todrop_cols).rename(columns=rename_cols)
    todrop_processed = updated_df["repeat_mapping"].isin(high_diff_repeats).index
    if solve_strat == "drop":
        updated_df = updated_df.drop(index=todrop_processed)
    # drop the repeats and concatenate with the filtered & updated values
    df = pd.concat(
        [
            df.drop(index=repeat_subset.index).assign(might_rancemic=lambda x: [False] * len(x)),
            updated_df.assign(might_rancemic=lambda x: [True] * len(x)),
        ],
        ignore_index=True,
    )
    logger.info(f"Final number of points: {len(df)}")
    return df, final_cols
