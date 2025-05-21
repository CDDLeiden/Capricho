"""Collection of functions to flag compounds based on specific criteria.

Not all the functions are annotating compounds to be removed from the dataset. Some are
used to annotate processing steps that occurred during the data processing pipeline, like:

- Salt/solvent removal (applied to the canonical SMILES since they're kept as-is by CompoundMapper)
- Calculated pChEMBL value (used when the this measure is absent in ChEMBL and calculated from nM ... etc readouts)
- Potential duplicates (when this one is found, CompoundMapper will keep only one of those in the dataset. If
  the user wants to investigate, please pass the `keep_duplicates` flag to True)

"""

import pandas as pd

from ..core.pandas_helper import add_comment


### Marking readouts that will be DROPPED by CompoundMapper ###
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


def flag_to_remove_mixture_compounds(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'mixture_compounds' is True."""
    return add_comment(
        df,
        comment="Mixture in SMILES",
        criteria_func=lambda x: x.str.contains(".", regex=False),
        target_column="standard_smiles",
        comment_type="d",
    )


def flag_undefined_stereochemistry(df: pd.DataFrame) -> pd.DataFrame:
    """Mark compounds with undefined stereochemistry based on a predefined boolean mask."""
    return add_comment(
        df,
        comment="Undefined stereochemistry",
        criteria_func=lambda x: x > 0,
        target_column="undefined_stereocenters",
        comment_type="d",
    )


### Marking readouts that were PROCESSED by CompoundMapper ###
def flag_calculated_pchembl(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows where 'calculated_pchembl' is True."""
    return add_comment(
        df,
        comment="Calculated pChEMBL",
        comment_type="p",
    )


def flag_salt_or_solvent_removal(df: pd.DataFrame) -> pd.DataFrame:
    """Marks rows with a salt/mixture on the canonical SMILES (not modified by CompoundMapper)"""
    return add_comment(
        df,
        comment="Salt/solvent removed",
        criteria_func=lambda x: x.str.contains(".", regex=False),
        target_column="canonical_smiles",
        comment_type="p",
    )


def flag_duplication(df: pd.DataFrame, dupli_id_subset: list[str]) -> pd.DataFrame:
    """Marks rows with a potential duplication after SMILES standardization & salt removal.

    Args:
        df: DataFrame to be processed.
        dupli_id_subset: subset of columns to be used for identifying duplicates.

    Returns:
        pd.DataFrame: DataFrame with duplicates marked in `data_processing_comment`
                      (or `data_dropping_comment` if comment_type='d' was used).
    """
    return (
        df.assign(temp_dupli_flag=lambda x: x.duplicated(subset=dupli_id_subset, keep=False))
        .pipe(
            add_comment,
            comment="Potential Duplication - same pChEMBL value reported for the same compound",
            criteria_func=lambda x: x,
            target_column="temp_dupli_flag",
            comment_type="p",
        )
        .drop(columns=["temp_dupli_flag"])
    )
