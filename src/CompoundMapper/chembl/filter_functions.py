import pandas as pd

from ..core.pandas_helper import add_comment


def flag_missing_canonical_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'canonical_smiles' is missing."""
    return add_comment(
        df,
        comment="Missing SMILES",
        criteria_func=lambda x: x.isna(),  # x is the Series df["canonical_smiles"]
        target_column="canonical_smiles",
        comment_type="d",
    )


def flag_missing_standard_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'standard_smiles' is missing."""
    return add_comment(
        df,
        comment="Missing Standard SMILES",
        criteria_func=lambda x: x.isna(),
        target_column="standard_smiles",
        comment_type="d",
    )


def flag_with_data_validity_comment(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'data_validity_comment' is present (not NA)."""
    return add_comment(
        df,
        comment="Data Validity Comment Present",
        criteria_func=lambda x: ~x.isna(),
        target_column="data_validity_comment",
        comment_type="d",
    )


def flag_potential_duplicate(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'potential_duplicate' is 0."""
    return add_comment(
        df,
        comment="Potential Duplicate",
        criteria_func=lambda x: x == 1,
        target_column="potential_duplicate",
        comment_type="d",
    )


def flag_calculated_pchembl(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'calculated_pchembl' is True."""
    return add_comment(
        df,
        comment="Calculated pChEMBL",
        comment_type="p",
    )


def rm_by_index(df, index, comment):
    """
    Flexible function to remove rows based on a list of pd.DataFrame indices.

    Args:
        df (pd.DataFrame): The dataframe to remove rows from.
        index (list): The list of indices to remove.
        comment (str): The comment to add to the removed rows.

    Returns:
        pd.DataFrame: The modified dataframe with the specified rows removed.
    """
    df = add_comment(
        df,
        comment=comment,
        index=index,
        criteria_func=lambda x: x.index.isin(index),
    )
    return df
