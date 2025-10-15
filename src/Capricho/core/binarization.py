"""Contain functions for binarizing bioactivity data; handling censored data and validating agreement between discrete and censored measurements"""

import numpy as np
import pandas as pd

from ..core.pandas_helper import add_comment
from ..logger import logger


def _invert_relation_for_pchembl(relation: str) -> str:
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
    censored_relation: str,
    threshold: float,
) -> bool:
    """Check if discrete and censored measurements agree for the same compound-target pair.

    Args:
        discrete_value: pchembl_value from "=" or "~" measurement
        censored_value: pchembl_value from censored measurement
        censored_relation: The censored relation ("<", ">", "<<", ">>", "<=", ">=")
        threshold: Activity threshold for binarization

    Returns:
        True if measurements agree, False otherwise
    """
    # Invert the censored relation for pchembl comparison
    pchembl_relation = _invert_relation_for_pchembl(censored_relation)

    # Check if censored measurement is consistent with discrete measurement
    if pchembl_relation in [">", ">>"]:
        # Censored says "pchembl > censored_value" (active)
        # Discrete should have pchembl >= censored_value
        return discrete_value >= censored_value
    elif pchembl_relation in ["<", "<<"]:
        # Censored says "pchembl < censored_value" (inactive)
        # Discrete should have pchembl <= censored_value
        return discrete_value <= censored_value
    else:
        # Should not reach here for censored relations
        return True


def binarize_aggregated_data(
    df: pd.DataFrame,
    threshold: float = 6.0,
    value_column: str = "pchembl_value_mean",
    compound_id_col: str = "connectivity",
    target_id_col: str = "target_chembl_id",
    relation_col: str = "standard_relation",
    output_binary_col: str = "activity_binary",
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

    Returns:
        DataFrame with binary activity labels and conflict flags
    """
    required_cols = [compound_id_col, target_id_col, value_column]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    df = df.copy()

    # Check if standard_relation column exists
    has_relation_col = relation_col in df.columns
    if not has_relation_col:
        logger.warning(f"Column '{relation_col}' not found. Assuming all measurements have '=' relation.")
        df[relation_col] = "="

    # Create mapping for same compound-target pairs (ignoring standard_relation)
    groupby_cols = [compound_id_col, target_id_col]

    # Initialize output columns
    df[output_binary_col] = np.nan
    conflict_flags = []

    # Process each compound-target group
    for group_key, group_df in df.groupby(groupby_cols):
        group_indices = group_df.index.tolist()
        relations_in_group = group_df[relation_col].unique()

        # Check for mixed relations
        has_discrete = "=" in relations_in_group or "~" in relations_in_group
        has_censored = any(rel in ["<", ">", "<=", ">=", ">>", "<<"] for rel in relations_in_group)

        if has_discrete and has_censored:
            # Need to validate agreement between discrete and censored measurements
            discrete_rows = group_df[group_df[relation_col].isin(["=", "~"])]
            censored_rows = group_df[group_df[relation_col].isin(["<", ">", "<=", ">=", ">>", "<<"])]

            # Check all discrete-censored pairs for agreement
            agreements = []
            for _, discrete_row in discrete_rows.iterrows():
                discrete_val = discrete_row[value_column]
                if pd.isna(discrete_val):
                    continue

                for _, censored_row in censored_rows.iterrows():
                    censored_val = censored_row[value_column]
                    censored_rel = censored_row[relation_col]
                    if pd.isna(censored_val):
                        continue

                    agrees = _check_measurement_agreement(discrete_val, censored_val, censored_rel, threshold)
                    agreements.append(agrees)

            # Flag if any disagreements found
            if agreements and not all(agreements):
                conflict_flags.extend(group_indices)
                logger.debug(
                    f"Conflict detected for {compound_id_col}={group_key[0]}, "
                    f"{target_id_col}={group_key[1]}: discrete and censored measurements disagree"
                )

        # Apply binarization logic to each row in the group
        for idx, row in group_df.iterrows():
            value = row[value_column]
            relation = row[relation_col]

            if pd.isna(value):
                df.loc[idx, output_binary_col] = np.nan
                continue

            # Apply binarization based on relation type
            if relation == "=":
                # Discrete measurement: compare directly
                df.loc[idx, output_binary_col] = 1 if value >= threshold else 0

            elif relation == "~":
                # Approximate measurement (±0.5 log units uncertainty)
                # Use the lower bound for conservative classification
                lower_bound = value - 0.5
                df.loc[idx, output_binary_col] = 1 if lower_bound >= threshold else 0

            elif relation in ["<", "<=", "<<"]:
                # Censored active: compound is MORE active than reported value
                # Since pchembl is -log, higher value = more active
                # If reported pchembl >= threshold, definitely active
                df.loc[idx, output_binary_col] = 1 if value >= threshold else 0

            elif relation in [">", ">=", ">>"]:
                # Censored inactive: compound is LESS active than reported value
                # Since pchembl is -log, lower value = less active
                # If reported pchembl <= threshold, definitely inactive
                df.loc[idx, output_binary_col] = 0 if value <= threshold else 1

            else:
                logger.warning(f"Unknown relation '{relation}' at index {idx}. Skipping binarization.")
                df.loc[idx, output_binary_col] = np.nan

    # Flag rows with conflicts
    if conflict_flags:
        logger.info(
            f"Found {len(conflict_flags)} measurements with disagreements between "
            "discrete (=) and censored (<, >) values"
        )
        df = add_comment(
            df=df,
            comment="Non-agreeing discrete and censored values",
            criteria_func=lambda x: x.index.isin(conflict_flags),
            target_column=value_column,
            comment_type="d",
        )

    # Convert binary column to Int64 (nullable integer)
    df[output_binary_col] = df[output_binary_col].astype("Int64")

    logger.info(
        f"Binarization complete. Threshold: {threshold}. "
        f"Active (1): {(df[output_binary_col] == 1).sum()}, "
        f"Inactive (0): {(df[output_binary_col] == 0).sum()}, "
        f"Missing: {df[output_binary_col].isna().sum()}"
    )

    return df
