"""Module with method to get similar compounds from the ChEMBL api using SMILES strings."""

from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from chemFilters.chem.standardizers import ChemStandardizer

from ..core.fp_utils import calculate_mixed_FPs
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import repeated_indices_from_array_series
from ..logger import logger


def get_and_curate_chembl_compounds(
    smiles: list[str], similarity: float, n_threads: int = 5, chirality: bool = True
) -> pd.DataFrame:
    """Get similar compounds from the ChEMBL API using multiple threads and identify
    repeats based on fingerprints. API call is forced to a limit of 5 calls per second
    not to overload ChEMBL servers. The fingerprints used to identify repeats is a combination of
    rdkit and morgan. Compounds identified as repeats will have the smallest `molecule_chembl_id`
    from the identified repeats listed in the `repeats` column.

    Returned fields (other than the ChEMBL fields):
    - querySmiles -> the original SMILES used to query the similar compounds
    - standard_smiles -> the standardized SMILES of the compound
    - repeats -> the smallest `molecule_chembl_id` from the identified repeats (based on fingerprint
        similarity after compound standardization + salt/solvent removal)

    Usage:
    >>> from CompoundMapper.chembl.similarity import get_and_curate_chembl_compounds
    >>> aspirin = 'O=C(C)Oc1ccccc1C(=O)O'
    >>> df = get_and_curate_chembl_compounds([aspirin], similarity=70)

    Args:
        smiles: list of SMILES to find similar molecules to.
        similarity: similarity threshold to use for the search. Value should be between 40 and 100.
        n_threads: Number of threads to use for searching the similar compounds. Defaults to 5.
        chirality: Whether to consider chirality when identifying repeats. Defaults to True.

    Returns:
        pd.DataFrame: a DataFrame with the processed similar molecules.
    """
    from .api.webresource import get_similarity_compound_table

    if not isinstance(smiles, list):
        raise TypeError("The 'smiles' argument must be a list of SMILES.")

    extracted = []

    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for smi in smiles:
            futures.append(executor.submit(get_similarity_compound_table, smi, similarity))
        for future in as_completed(futures):
            extracted.append(future.result())

    similar_compounds = pd.concat(extracted, ignore_index=True)

    if similar_compounds.empty:
        logger.warning(
            "No similar compounds found for the provided SMILES with the given similarity threshold {similarity}."
        )
        return pd.DataFrame()

    stdzer = ChemStandardizer(from_smi=True, n_jobs=8, verbose=False)
    df = (
        similar_compounds.assign(standard_smiles=lambda x: stdzer(x["canonical_smiles"]))
        .dropna(subset=["standard_smiles"])  # drop if no structure is found
        .query("standard_smiles.notna()")
        .assign(final_smiles=lambda x: x["standard_smiles"].apply(clean_mixtures))
        .drop(columns="standard_smiles")
        .rename(columns={"final_smiles": "standard_smiles"})
    )

    if df.shape[0] == 1:
        return df
    elif df.shape[0] == 0:
        return pd.DataFrame()

    # Calculate fingerprints and identify repeats
    fps = calculate_mixed_FPs(
        df["standard_smiles"].tolist(), n_jobs=n_threads, morgan_kwargs={"useChirality": chirality}
    )
    df = df.assign(fps=fps).assign(repeats=False)
    repeats_idxs = repeated_indices_from_array_series(df["fps"])
    for repeats in repeats_idxs:
        min_id = sorted(df.loc[repeats, "molecule_chembl_id"], key=lambda _id: int(_id.lstrip("CHEMBL")))[0]
        df.loc[repeats, "repeats"] = min_id

    df = df.drop(columns=["fps"])

    return df
