"""Module containing helper functions for manipulating pandas DataFrames"""

import functools
from typing import Union

import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation

from ..logger import logger


def format_value(x) -> str:
    """Helper function to format a value to a string with 4 decimal places. This function
    is used to store the original pChEMBL values as strings separated by ";".

    Args:
        x: Value to be formatted

    Returns:
        x: Formatted value
    """
    if isinstance(x, float):
        if x.is_integer() and 1000 <= x <= 9999:
            return f"{int(x)}"  # Return as integer string if it's a year
        return f"{x:.2f}"
    elif isinstance(x, int):
        return str(x)
    else:
        return x


def aggr_val_series(series: pd.Series) -> str:
    """ "Aggregate a pandas Series into a string with values separated by a semicolon."""
    return ";".join([format_value(x) for x in series])


def get_mad(values) -> Union[float, np.float64]:
    """Calculate the MAD for a list of numerical values. If only one value, return NaN."""
    if len(values) > 1:
        return median_abs_deviation(values)
    else:
        return np.nan


def merge_dataframes(dfs, id_cols) -> pd.DataFrame:
    """
    Merge a list of DataFrames based on id_cols. Useful reference for merges:
    https://stackoverflow.com/questions/53645882/pandas-merging-101

    Args:
        dfs: List of DataFrames to be merged
        id_cols: List of columns to be used as keys for merging

    Returns:
        merged_df: Merged DataFrame
    """
    return functools.reduce(lambda left, right: pd.merge(left, right, on=id_cols, how="inner"), dfs)


def apply_func_grpd(grpd, func: callable, idcols: list, *cols: list) -> pd.DataFrame:
    """Apply a function to a list of columns (*cols) of a grouped dataframe and
    merge the results based on id_cols.

    Args:
        grpd: grouped dataframe
        func: callable function to be applied
        idcols: list of columns to be used as index
    """
    results = []
    for col in cols:
        try:
            results.append(grpd[col].apply(func, group_keys=True).reset_index().set_index(idcols).copy())
        except Exception as e:
            logger.error(f"Error applying function to column {col}: {e}")
            logger.error(f"Exception type: {type(e).__name__}")
            raise e
    return pd.concat(results, ignore_index=False, axis=1).reset_index()


def assign_stats(df: pd.DataFrame, sep=";", value_col="pchembl_value") -> pd.DataFrame:
    """Assign statistics to a DataFrame based on a column with multiple values separated by
    a particular separator, e.g. ';'.

    Args:
        df: pd.DataFrame to be processed.
        sep: string separating the values. Defaults to ';'.
        value_col: column containing the values to be processed. Defaults to "pchembl_value".

    Returns:
        pd.DataFrame: DataFrame with the statistics assigned. as columns:
    >>> new_cols = [
    >>>     f"{value_col}_mean",
    >>>     f"{value_col}_std",
    >>>     f"{value_col}_median",
    >>>     f"{value_col}_mad",
    >>>     f"{value_col}_counts",
    >>> ]
    """

    value_series = df[value_col].astype(str).str.split(sep).apply(lambda x: list(map(float, x)))
    new_cols = [f"{value_col}{suffix}" for suffix in ["_mean", "_std", "_median", "_mad", "_counts"]]
    values = [
        value_series.apply(np.mean),
        value_series.apply(np.std),
        value_series.apply(np.median),
        value_series.apply(get_mad),
        value_series.apply(lambda x: len(x)),
    ]
    for c, v in zip(new_cols, values):
        df[c] = v
    return df


def find_dict_in_dataframe(df):
    cols_w_dicts = []
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, dict)).any():
            logger.info(f"Column '{col}' contains dictionaries.")
            logger.info(
                "Rows with dictionaries: "
                + " ".join(df[df[col].apply(lambda x: isinstance(x, dict))].index.astype(str))
            )
            cols_w_dicts.append(col)
    if cols_w_dicts:
        return cols_w_dicts
