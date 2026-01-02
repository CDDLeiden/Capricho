"""Module for cleaning up specified data quality flags, and transforming the bioactivity data
into matrices where rows are compounds and columns are affinity values to be used for multitask modeling.
"""

from typing import List, Optional

import pandas as pd

from ..core.default_fields import DATA_DROPPING_COMMENT
from ..logger import logger


def prepare_multitask_data(
    df: pd.DataFrame,
    task_col: str,
    value_col: str,
    compound_col: str,
    smiles_col: str,
    remove_flags: Optional[List[str]] = None,
    id_columns: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Transform aggregated data to multitask format (activity matrix).

    This function pivots aggregated bioactivity data to create an activity matrix where:
    - Rows represent unique compounds (identified by compound_col)
    - Columns represent tasks (e.g., different targets)
    - Values are the bioactivity measurements (e.g., pchembl_value_mean)

    Args:
        df: Aggregated DataFrame from aggregate_data() with bioactivity statistics
        task_col: Column to use as task identifier (e.g., "target_chembl_id")
        value_col: Column containing values to pivot (e.g., "pchembl_value_mean")
        compound_col: Column for compound identity (e.g., "connectivity" or "smiles")
        smiles_col: Column containing SMILES strings
        remove_flags: List of quality flags to remove. Rows where data_dropping_comment
            contains any of these flags will be filtered out before pivoting.
            If None, no filtering is applied.
        id_columns: List of additional columns to combine with task_col for creating
            composite task identifiers. Use this when data was aggregated with
            --id-columns (e.g., ["assay_tissue"]) to prevent losing information.

    Returns:
        DataFrame with compounds as rows (indexed by compound_col), tasks as columns,
        and a smiles column. Missing values are represented as NaN.
    """
    df = df.copy()

    # Filter out flagged data if remove_flags is provided
    if remove_flags is not None and len(remove_flags) > 0:
        if DATA_DROPPING_COMMENT not in df.columns:
            logger.warning(
                f"Column '{DATA_DROPPING_COMMENT}' not found in DataFrame. " "No filtering will be applied."
            )
        else:
            initial_len = len(df)
            # Create a mask for rows containing any of the flags
            mask = pd.Series(False, index=df.index)
            for flag in remove_flags:
                flag_mask = df[DATA_DROPPING_COMMENT].fillna("").str.contains(flag, regex=False, na=False)
                mask = mask | flag_mask
                n_flagged = flag_mask.sum()
                if n_flagged > 0:
                    logger.info(f"Filtering {n_flagged} rows with flag: '{flag}'")

            df = df[~mask].copy()
            n_removed = initial_len - len(df)
            if n_removed > 0:
                logger.info(f"Total rows filtered out: {n_removed} ({n_removed/initial_len*100:.1f}%)")

    # Validate required columns exist
    required_cols = [compound_col, task_col, value_col, smiles_col]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    # Create composite task column if id_columns are provided
    effective_task_col = task_col
    if id_columns is not None and len(id_columns) > 0:
        missing_id_cols = [col for col in id_columns if col not in df.columns]
        if missing_id_cols:
            raise ValueError(f"Missing id_columns: {missing_id_cols}")

        # Create composite task identifier by joining task_col with id_columns
        composite_name = "_composite_task"
        df[composite_name] = df[task_col].astype(str)
        for col in id_columns:
            df[composite_name] = df[composite_name] + "-" + df[col].fillna("").astype(str)
        effective_task_col = composite_name
        logger.info(f"Created composite task column from: {task_col} + {id_columns}")

    # Check for duplicates (multiple values per compound-task pair)
    dup_check = df.groupby([compound_col, effective_task_col]).size()
    duplicates = dup_check[dup_check > 1]
    if len(duplicates) > 0:
        n_dup_pairs = len(duplicates)
        n_extra_rows = duplicates.sum() - n_dup_pairs
        logger.warning(
            f"Found {n_dup_pairs} compound-task pairs with multiple values ({n_extra_rows} extra rows). "
            f"Only the first value will be kept. "
            f"If your data was aggregated with --id-columns, use the same columns here via --id-columns."
        )

    # Pivot the data to create activity matrix
    logger.info(
        f"Creating activity matrix with {df[compound_col].nunique()} compounds "
        f"and {df[effective_task_col].nunique()} tasks"
    )

    activity_matrix = df.pivot_table(
        index=compound_col,
        columns=effective_task_col,
        values=value_col,
        aggfunc="first",
    )

    # Reset the columns name to remove the task_col label
    activity_matrix.columns.name = None

    # Add SMILES column back by taking the first SMILES for each compound
    smiles_map = df.groupby(compound_col)[smiles_col].first()
    activity_matrix[smiles_col] = smiles_map

    logger.info(
        f"Activity matrix shape: {activity_matrix.shape[0]} compounds x {activity_matrix.shape[1]} columns"
    )

    return activity_matrix
