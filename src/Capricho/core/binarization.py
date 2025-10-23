"""Contain functions for binarizing bioactivity data; handling censored data and validating agreement between discrete and censored measurements"""

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..core.pandas_helper import add_comment
from ..logger import logger


def _truncate_dataframe(df: pd.DataFrame, limit: int = 15) -> pd.DataFrame:
    """Truncate DataFrame values to a specified length for display purposes.

    Args:
        df: DataFrame to truncate
        limit: Maximum character length for each cell value

    Returns:
        DataFrame with truncated string values
    """
    if pd.__version__ > "2.1.0":  # applymap got deprecated in 2.1.0
        return df.map(lambda x: str(x)[:limit] + "..." if len(str(x)) > limit else str(x))
    else:
        return df.applymap(lambda x: str(x)[:limit] + "..." if len(str(x)) > limit else str(x))


def invert_relation_for_pchembl(relation: str) -> str:
    """Inverts comparison relation for pchembl values.

    Since pchembl = -log10(Molar), higher pchembl = more active (lower concentration).
    Therefore, standard_relation directions must be inverted:
    - standard_relation "<" (low concentration, active) → pchembl ">" (high value, active)
    - standard_relation ">" (high concentration, inactive) → pchembl "<" (low value, inactive)

    Args:
        relation: Original standard_relation from ChEMBL ("=", "<", ">", "<=", ">=", "~", ">>", "<<")

    Returns:
        Inverted relation for pchembl comparison
    """
    inversion_map = {
        "<": ">",
        ">": "<",
        "<=": ">=",
        ">=": "<=",
        "=": "=",
        "~": "~",
        ">>": "<<",
        "<<": ">>",
    }
    if relation not in inversion_map:
        raise ValueError(f"Unknown relation: {relation}")
    return inversion_map[relation]


def _check_measurement_agreement(
    discrete_value: float,
    censored_value: float,
    censored_pchembl_relation: str,
) -> bool:
    """Check if discrete and censored measurements agree for the same compound-target pair.

    Args:
        discrete_value: pchembl_value from "=" or "~" measurement
        censored_value: pchembl_value from censored measurement
        censored_pchembl_relation: The censored relation already converted to pchembl space
            (">", "<", ">>", "<<", ">=", "<=")

    Returns:
        True if measurements agree, False otherwise
    """
    # Check if censored measurement is consistent with discrete measurement
    if censored_pchembl_relation in [">", ">>"]:
        # Censored says "pchembl > censored_value" (active)
        # Discrete should have pchembl >= censored_value
        return discrete_value >= censored_value
    elif censored_pchembl_relation in ["<", "<<"]:
        # Censored says "pchembl < censored_value" (inactive)
        # Discrete should have pchembl <= censored_value
        return discrete_value <= censored_value
    else:
        # Should not reach here for censored relations
        return True


def _classify_by_relation(value: float, relation: str, threshold: float) -> int:
    """Classify a single measurement as active (1) or inactive (0) based on relation type.

    Args:
        value: pchembl_value to classify
        relation: standard_relation ("=", "~", "<", ">", "<=", ">=", "<<", ">>")
        threshold: Activity threshold for binarization

    Returns:
        1 for active, 0 for inactive
    """
    if relation == "=":
        return 1 if value >= threshold else 0

    elif relation == "~":
        lower_bound = value - 0.5
        return 1 if lower_bound >= threshold else 0

    elif relation in ["<", "<=", "<<"]:
        return 1 if value >= threshold else 0

    elif relation in [">", ">=", ">>"]:
        return 0 if value <= threshold else 1

    else:
        raise ValueError(f"Unknown relation: {relation}")


def _detect_conflicts(
    group_df: pd.DataFrame,
    value_column: str,
    pchembl_relation_col: str,
    output_binary_col: str,
) -> tuple[bool, str]:
    """Detect conflicts within a compound-target group.

    Args:
        group_df: DataFrame subset for one compound-target pair
        value_column: Column name with pchembl values
        pchembl_relation_col: Column name with pchembl_relation
        output_binary_col: Column name with binary labels

    Returns:
        Tuple of (has_conflict, conflict_type) where conflict_type is one of:
        "discrete_censored_disagreement", "binary_label_mismatch", or ""
    """
    relations_in_group = group_df[pchembl_relation_col].unique()
    has_discrete = "=" in relations_in_group or "~" in relations_in_group
    has_censored = any(rel in ["<", ">", "<=", ">=", ">>", "<<"] for rel in relations_in_group)

    if has_discrete and has_censored:
        discrete_rows = group_df[group_df[pchembl_relation_col].isin(["=", "~"])]
        censored_rows = group_df[group_df[pchembl_relation_col].isin(["<", ">", "<=", ">=", ">>", "<<"])]

        agreements = []
        for _, discrete_row in discrete_rows.iterrows():
            discrete_val = discrete_row[value_column]
            if pd.isna(discrete_val):
                continue

            for _, censored_row in censored_rows.iterrows():
                censored_val = censored_row[value_column]
                censored_pchembl_rel = censored_row[pchembl_relation_col]
                if pd.isna(censored_val):
                    continue

                agrees = _check_measurement_agreement(discrete_val, censored_val, censored_pchembl_rel)
                agreements.append(agrees)

        if agreements and not all(agreements):
            return True, "discrete_censored_disagreement"

    binary_labels = group_df[output_binary_col].dropna()
    if len(binary_labels) > 1 and len(binary_labels.unique()) > 1:
        return True, "binary_label_mismatch"

    return False, ""


def _generate_conflict_details(
    df: pd.DataFrame,
    conflict_indices: list,
    compound_id_col: str,
    target_id_col: str,
    pchembl_relation_col: str,
    value_column: str,
    output_binary_col: str,
    threshold: float,
) -> list[dict]:
    """Generate detailed conflict information for each compound-target pair.

    Args:
        df: DataFrame with binarized data
        conflict_indices: List of row indices with conflicts
        compound_id_col: Column name for compound IDs
        target_id_col: Column name for target IDs
        pchembl_relation_col: Column name for pchembl_relation
        value_column: Column name for pchembl values
        output_binary_col: Column name for binary labels
        threshold: Binarization threshold

    Returns:
        List of conflict detail dictionaries
    """
    if not conflict_indices:
        return []

    conflict_subset = df.loc[conflict_indices].copy()
    groupby_cols = [compound_id_col, target_id_col]
    if "mutation" in conflict_subset.columns:
        groupby_cols.append("mutation")

    conflict_details = []
    for group_key, group_df in conflict_subset.groupby(groupby_cols):
        has_conflict, conflict_type = _detect_conflicts(
            group_df, value_column, pchembl_relation_col, output_binary_col
        )

        if not has_conflict:
            continue

        measurements = []
        for _, row in group_df.iterrows():
            measurement = {
                "value": float(row[value_column]) if not pd.isna(row[value_column]) else None,
                "pchembl_relation": row[pchembl_relation_col],
                "binary": int(row[output_binary_col]) if not pd.isna(row[output_binary_col]) else None,
                "threshold_distance": (
                    float(row[value_column] - threshold) if not pd.isna(row[value_column]) else None
                ),
            }

            if "standard_relation" in row and not pd.isna(row["standard_relation"]):
                measurement["standard_relation"] = row["standard_relation"]
            if "assay_chembl_id" in row:
                measurement["assay"] = row["assay_chembl_id"]
            if "molecule_chembl_id" in row:
                measurement["molecule"] = row["molecule_chembl_id"]

            measurements.append(measurement)

        conflict_detail = {
            "compound_id": group_key[0],
            "target_id": group_key[1],
            "conflict_type": conflict_type,
            "measurements": measurements,
            "threshold": threshold,
        }

        if len(groupby_cols) > 2:
            conflict_detail["mutation"] = group_key[2]

        discrete_measurements = [
            m for m in measurements if m["pchembl_relation"] in ["=", "~"] and m["value"] is not None
        ]
        censored_measurements = [
            m
            for m in measurements
            if m["pchembl_relation"] in ["<", ">", "<=", ">=", "<<", ">>"] and m["value"] is not None
        ]

        if discrete_measurements and censored_measurements:
            discrete_parts = [f"{m['pchembl_relation']} {m['value']:.2f}" for m in discrete_measurements]
            censored_parts = [f"{m['pchembl_relation']} {m['value']:.2f}" for m in censored_measurements]

            discrete_str = ", ".join(discrete_parts)
            censored_str = ", ".join(censored_parts)

            conflict_detail["explanation"] = (
                f"Discrete measurement(s) ({discrete_str}) vs censored ({censored_str}). "
                f"Threshold: {threshold}"
            )
        else:
            binary_labels = [m["binary"] for m in measurements if m["binary"] is not None]
            conflict_detail["explanation"] = (
                f"Multiple measurements with different binary labels: {set(binary_labels)}"
            )

        conflict_details.append(conflict_detail)

    return conflict_details


def _log_and_flag_conflicts(
    df: pd.DataFrame,
    conflict_indices: list,
    compound_id_col: str,
    target_id_col: str,
    relation_col: str,
    value_column: str,
    output_binary_col: str,
) -> pd.DataFrame:
    """Log conflict details and flag conflicting rows in the DataFrame.

    Args:
        df: DataFrame with binarized data
        conflict_indices: List of row indices with conflicts
        compound_id_col: Column name for compound IDs
        target_id_col: Column name for target IDs
        relation_col: Column name for relations
        value_column: Column name for pchembl values
        output_binary_col: Column name for binary labels

    Returns:
        DataFrame with conflicts flagged via add_comment()
    """
    if not conflict_indices:
        return df

    logger.warning(
        f"Found {len(conflict_indices)} measurements with disagreements. "
        "These compound-target pairs have inconsistent measurements across different relation types."
    )

    conflict_subset = df.loc[conflict_indices].copy()

    logging_cols = [compound_id_col, target_id_col, relation_col, value_column, "mutation"]
    optional_cols = ["molecule_chembl_id", "assay_chembl_id", output_binary_col, "standard_relation"]
    for col in optional_cols:
        if col in conflict_subset.columns:
            logging_cols.append(col)

    logging_cols = [col for col in logging_cols if col in conflict_subset.columns]
    conflict_display = conflict_subset[logging_cols].sort_values(
        by=[target_id_col, compound_id_col, value_column], ascending=[True, True, False]
    )

    truncated_df = _truncate_dataframe(conflict_display, limit=15)
    logger.warning(
        f"Sample of conflicting measurements (showing up to 20 rows):\n{truncated_df.head(20).to_string(index=False)}"
    )

    df = add_comment(
        df=df,
        comment="Non-agreeing discrete and censored values",
        criteria_func=lambda x: x.index.isin(conflict_indices),
        target_column=value_column,
        comment_type="d",
    )

    return df


def save_conflict_report(
    conflict_details: list[dict],
    output_path: str | Path,
    threshold: float,
) -> None:
    """Save conflict report to JSON file.

    Args:
        conflict_details: List of conflict detail dictionaries
        output_path: Path to save the JSON file
        threshold: Binarization threshold used
    """
    report = {
        "summary": {
            "total_conflicts": len(conflict_details),
            "discrete_censored_disagreements": sum(
                1 for c in conflict_details if c["conflict_type"] == "discrete_censored_disagreement"
            ),
            "binary_label_mismatches": sum(
                1 for c in conflict_details if c["conflict_type"] == "binary_label_mismatch"
            ),
            "threshold": threshold,
        },
        "conflicts": conflict_details,
    }

    output_path = Path(output_path)
    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)

    logger.info(f"Conflict report saved to {output_path}")


def binarize_aggregated_data(
    df: pd.DataFrame,
    threshold: float = 6.0,
    value_column: str = "pchembl_value_mean",
    compound_id_col: str = "connectivity",
    target_id_col: str = "target_chembl_id",
    relation_col: str = "standard_relation",
    output_binary_col: str = "activity_binary",
    compare_across_mutants: bool = False,
    conflict_report_path: Optional[str | Path] = None,
) -> pd.DataFrame:
    """Binarize aggregated bioactivity data based on activity threshold and standard_relation.

    This function converts continuous pchembl values to binary labels (0=inactive, 1=active)
    while properly handling censored measurements and approximate values, and validating
    agreement between discrete and censored measurements for the same compound-target pair.

    Key logic:
    - standard_relation "=": compare value to threshold directly
    - standard_relation "~": approximate (±0.5 log units); uses lower bound for conservative classification
    - standard_relation "<", "<<" (low concentration): if pchembl >= threshold → active (1)
    - standard_relation ">", ">>" (high concentration): if pchembl <= threshold → inactive (0)
    - Mixed relations: validate agreement and flag conflicts

    Args:
        df: Aggregated DataFrame from aggregate_data() with pchembl statistics
        threshold: Activity threshold for binarization (default 6.0 = 1 µM)
        value_column: Which aggregated column to use (default: "pchembl_value_mean")
        compound_id_col: Column identifying compounds (default: "connectivity")
        target_id_col: Column identifying targets (default: "target_chembl_id")
        relation_col: Column with standard_relation values (default: "standard_relation")
        output_binary_col: Name for output binary column (default: "activity_binary")
        compare_across_mutants: If False (default), different mutations are treated as separate
            compound-target pairs for conflict detection. If True, measurements on different
            mutants are compared and flagged if they disagree.
        conflict_report_path: Optional path to save detailed conflict report as JSON

    Returns:
        DataFrame with binary activity labels, pchembl_relation column, and conflict flags
    """
    required_cols = [compound_id_col, target_id_col, value_column]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    df = df.copy()

    if relation_col not in df.columns:
        logger.warning(f"Column '{relation_col}' not found. Assuming all measurements have '=' relation.")
        df[relation_col] = "="

    pchembl_relation_col = "pchembl_relation"
    df[pchembl_relation_col] = df[relation_col].apply(invert_relation_for_pchembl)

    groupby_cols = [compound_id_col, target_id_col]
    if not compare_across_mutants and "mutation" in df.columns:
        groupby_cols.append("mutation")

    df[output_binary_col] = np.nan
    conflict_indices = []

    for idx, row in df.iterrows():
        value = row[value_column]
        relation = row[relation_col]

        if pd.isna(value):
            continue

        try:
            df.loc[idx, output_binary_col] = _classify_by_relation(value, relation, threshold)
        except ValueError as e:
            logger.warning(f"{e} at index {idx}. Skipping binarization.")

    for group_key, group_df in df.groupby(groupby_cols):
        has_conflict, conflict_type = _detect_conflicts(
            group_df, value_column, pchembl_relation_col, output_binary_col
        )

        if has_conflict:
            conflict_indices.extend(group_df.index.tolist())
            logger.debug(
                f"Conflict ({conflict_type}) for {compound_id_col}={group_key[0]}, "
                f"{target_id_col}={group_key[1]}"
            )

    df = _log_and_flag_conflicts(
        df,
        conflict_indices,
        compound_id_col,
        target_id_col,
        pchembl_relation_col,
        value_column,
        output_binary_col,
    )

    if conflict_report_path and conflict_indices:
        conflict_details = _generate_conflict_details(
            df,
            conflict_indices,
            compound_id_col,
            target_id_col,
            pchembl_relation_col,
            value_column,
            output_binary_col,
            threshold,
        )
        save_conflict_report(conflict_details, conflict_report_path, threshold)

    df[output_binary_col] = df[output_binary_col].astype("Int64")

    logger.info(
        f"Binarization complete. Threshold: {threshold}. "
        f"Active (1): {(df[output_binary_col] == 1).sum()}, "
        f"Inactive (0): {(df[output_binary_col] == 0).sum()}, "
        f"Missing: {df[output_binary_col].isna().sum()}"
    )

    return df
