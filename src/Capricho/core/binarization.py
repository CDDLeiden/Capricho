"""Contain functions for binarizing bioactivity data; handling censored data and validating agreement between discrete and censored measurements"""

import json
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..core.pandas_helper import add_comment
from ..logger import logger


def _calculate_mcc(tp: int, tn: int, fp: int, fn: int) -> float:
    """Calculate Matthews Correlation Coefficient from confusion matrix values.

    MCC = (TP*TN - FP*FN) / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))

    Args:
        tp: True positives
        tn: True negatives
        fp: False positives
        fn: False negatives

    Returns:
        MCC value between -1.0 and 1.0, or 0.0 if denominator is zero
    """
    numerator = (tp * tn) - (fp * fn)
    denominator_parts = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)

    if denominator_parts == 0:
        return 0.0

    import math

    denominator = math.sqrt(denominator_parts)

    return numerator / denominator


def _calculate_assay_compatibility_mcc(
    df: pd.DataFrame,
    compound_id_col: str,
    target_id_col: str,
    output_binary_col: str,
) -> float:
    """Calculate MCC measuring agreement between measurements for same compound-target pairs.

    Uses pairwise comparisons: for each compound-target pair with multiple measurements,
    compares all pairs of binary labels to build confusion matrix, then calculates MCC.

    Args:
        df: DataFrame with binarized data
        compound_id_col: Column name for compound IDs
        target_id_col: Column name for target IDs
        output_binary_col: Column name with binary labels

    Returns:
        MCC value between -1.0 and 1.0
    """
    groupby_cols = [compound_id_col, target_id_col]
    if "mutation" in df.columns:
        groupby_cols.append("mutation")

    tp = tn = fp = fn = 0

    for _, group_df in df.groupby(groupby_cols):
        binary_labels = group_df[output_binary_col].dropna().tolist()

        # Skip if less than 2 measurements
        if len(binary_labels) < 2:
            continue

        # Pairwise comparisons within this group
        for i in range(len(binary_labels)):
            for j in range(i + 1, len(binary_labels)):
                label_i = binary_labels[i]
                label_j = binary_labels[j]

                # Treat label_i as "truth" and label_j as "prediction"
                if label_i == 1 and label_j == 1:
                    tp += 1
                elif label_i == 0 and label_j == 0:
                    tn += 1
                elif label_i == 0 and label_j == 1:
                    fp += 1
                elif label_i == 1 and label_j == 0:
                    fn += 1

    return _calculate_mcc(tp, tn, fp, fn)


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
    output_binary_col: str,
) -> tuple[bool, str]:
    """Detect conflicts within a compound-target group based on binarization outcomes.

    A conflict occurs when multiple measurements for the same compound-target pair
    result in different binary labels (active vs inactive).

    Args:
        group_df: DataFrame subset for one compound-target pair
        output_binary_col: Column name with binary labels

    Returns:
        Tuple of (has_conflict, conflict_type) where conflict_type is "binary_label_mismatch" or ""
    """
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
        has_conflict, conflict_type = _detect_conflicts(group_df, output_binary_col)

        if not has_conflict:
            continue

        measurements = []
        for _, row in group_df.iterrows():
            measurement = {
                "value": float(row[value_column]) if not pd.isna(row[value_column]) else None,
                "pchembl_relation": row[pchembl_relation_col],
                "binary": int(row[output_binary_col]) if not pd.isna(row[output_binary_col]) else None,
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

        # Group measurements by binary outcome
        active_measurements = [m for m in measurements if m["binary"] == 1]
        inactive_measurements = [m for m in measurements if m["binary"] == 0]

        # Create vote summary
        vote_summary = {
            "active_votes": len(active_measurements),
            "inactive_votes": len(inactive_measurements),
        }
        conflict_detail["vote_summary"] = vote_summary

        # Build explanation with vote counts and measurement details
        explanation_parts = []

        if active_measurements:
            active_strs = [f"{m['pchembl_relation']}{m['value']:.2f}" for m in active_measurements]
            explanation_parts.append(f"Active ({len(active_measurements)} votes): {', '.join(active_strs)}")

        if inactive_measurements:
            inactive_strs = [f"{m['pchembl_relation']}{m['value']:.2f}" for m in inactive_measurements]
            explanation_parts.append(
                f"Inactive ({len(inactive_measurements)} votes): {', '.join(inactive_strs)}"
            )

        conflict_detail["explanation"] = " | ".join(explanation_parts) + f" | Threshold: {threshold}"

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
        has_conflict, conflict_type = _detect_conflicts(group_df, output_binary_col)

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

    # Calculate summary statistics
    n_active = (df[output_binary_col] == 1).sum()
    n_inactive = (df[output_binary_col] == 0).sum()
    n_missing = df[output_binary_col].isna().sum()
    n_total = len(df)

    # Calculate conflict statistics
    n_conflicting_measurements = len(conflict_indices)

    # Count unique compound-target pairs with conflicts
    if conflict_indices:
        conflict_subset = df.loc[conflict_indices]
        conflict_groupby_cols = [compound_id_col, target_id_col]
        if not compare_across_mutants and "mutation" in conflict_subset.columns:
            conflict_groupby_cols.append("mutation")
        n_conflicting_pairs = conflict_subset.groupby(conflict_groupby_cols).ngroups
    else:
        n_conflicting_pairs = 0

    # Calculate MCC for assay compatibility (only for pairs with multiple measurements)
    mcc = _calculate_assay_compatibility_mcc(df, compound_id_col, target_id_col, output_binary_col)

    # Log comprehensive summary
    logger.info(
        f"BINARIZATION SUMMARY\n"
        f"Threshold: {threshold}\n"
        f"\nBinary labels:\n"
        f"  Active (1):   {n_active:>6} ({n_active/n_total*100:>5.1f}%)\n"
        f"  Inactive (0): {n_inactive:>6} ({n_inactive/n_total*100:>5.1f}%)\n"
        f"  Missing:      {n_missing:>6} ({n_missing/n_total*100:>5.1f}%)\n"
        f"\nConflict analysis:\n"
        f"  Measurements with conflicts:      {n_conflicting_measurements:>6}\n"
        f"  Compound-target pairs affected:   {n_conflicting_pairs:>6}\n"
        f"  Datapoints lost if conflicts removed: {n_conflicting_measurements:>6} "
        f"({n_conflicting_measurements/n_total*100:>5.1f}%)\n"
        f"\nAssay compatibility:\n"
        f"  MCC (agreement between measurements): {mcc:>6.3f}\n"
    )

    return df
