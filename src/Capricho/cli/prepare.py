"""Clean aggregated bioactivity data by removing quality flags and duplicates,
and pivot into activity matrices for multitask modeling.
"""

from typing import List, Optional

import pandas as pd

from ..core.pandas_helper import filter_dropping_flags
from ..logger import logger


def clean_data(
    df: pd.DataFrame,
    drop_flags: Optional[List[str]] = None,
    deduplicate: bool = False,
    value_col: str = "pchembl_value",
) -> pd.DataFrame:
    """Clean aggregated bioactivity data by deduplicating and filtering quality flags.

    Orchestrates cleaning steps in the correct order:
    1. Deduplicate (if requested): removes identical values within aggregated rows
       and recalculates statistics (mean, median, std, counts).
    2. Drop flags: removes rows matching any of the specified quality flags.

    Appropriate flags for dropping include unit errors, undefined stereochemistry,
    assay size issues, and mixtures. For potential duplicates, prefer using
    ``deduplicate=True`` instead of dropping — this keeps one measurement while
    removing extras, rather than discarding the entire row.

    Annotation error resolution is a separate step handled by
    ``resolve_annotation_errors()`` (available via ``capricho prepare
    --resolve-annotation-error``).

    Example::

        from Capricho.cli.prepare import clean_data, prepare_multitask_data

        cleaned = clean_data(df, drop_flags=["Unit Annotation Error"], deduplicate=True)
        matrix = prepare_multitask_data(cleaned, task_col="target_chembl_id", ...)

    Args:
        df: Aggregated DataFrame from aggregate_data().
        drop_flags: List of quality flags to remove. Rows where data_dropping_comment
            contains any of these flags will be filtered out.
        deduplicate: If True, remove duplicate values within aggregated rows
            and recalculate statistics.
        value_col: Column containing the activity values (e.g., "pchembl_value"
            or "standard_value"). Used for deduplication and stats recalculation.

    Returns:
        Cleaned DataFrame.
    """
    from ..analysis import deduplicate_aggregated_values, recalculate_aggregated_stats

    df = df.copy()

    # Step 1: Deduplicate
    if deduplicate:
        logger.info("Deduplicating identical values within aggregated rows...")
        initial_total = df[value_col].apply(
            lambda x: len(str(x).split("|")) if pd.notna(x) else 0
        ).sum()
        df = deduplicate_aggregated_values(df, value_column=value_col)
        final_total = df[value_col].apply(
            lambda x: len(str(x).split("|")) if pd.notna(x) else 0
        ).sum()
        logger.info(f"Deduplication removed {initial_total - final_total} duplicate values")

        logger.info("Recalculating statistics after deduplication...")
        df = recalculate_aggregated_stats(df, value_column=value_col)

    # Step 2: Drop flags
    if drop_flags:
        df = filter_dropping_flags(df, drop_flags)

    return df


def prepare_multitask_data(
    df: pd.DataFrame,
    task_col: str,
    value_col: str,
    compound_col: str,
    smiles_col: str,
    id_columns: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Transform aggregated data to multitask format (activity matrix).

    This function pivots aggregated bioactivity data to create an activity matrix where:
    - Rows represent unique compounds (identified by compound_col)
    - Columns represent tasks (e.g., different targets)
    - Values are the bioactivity measurements (e.g., pchembl_value_mean)

    Use ``clean_data()`` before calling this function to filter quality flags
    and deduplicate values.

    Args:
        df: Aggregated DataFrame from aggregate_data() with bioactivity statistics.
        task_col: Column to use as task identifier (e.g., "target_chembl_id").
        value_col: Column containing values to pivot (e.g., "pchembl_value_mean").
        compound_col: Column for compound identity (e.g., "connectivity" or "smiles").
        smiles_col: Column containing SMILES strings.
        id_columns: List of additional columns to combine with task_col for creating
            composite task identifiers. Use this when data was aggregated with
            --id-columns (e.g., ["assay_tissue"]) to prevent losing information.

    Returns:
        DataFrame with compounds as rows (indexed by compound_col), tasks as columns,
        and a smiles column. Missing values are represented as NaN.
    """
    df = df.copy()

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
    if smiles_col != compound_col:
        smiles_map = df.groupby(compound_col)[smiles_col].first()
        activity_matrix[smiles_col] = smiles_map

    logger.info(
        f"Activity matrix shape: {activity_matrix.shape[0]} compounds x {activity_matrix.shape[1]} columns"
    )

    return activity_matrix
