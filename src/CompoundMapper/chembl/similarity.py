"""Module with method to get similar compounds from the ChEMBL api using SMILES strings."""

from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

from .api import get_similarity_compound_table
from .rate_limit import rate_limit


@rate_limit(max_per_second=5)
def get_similars_from_smiles(smiles: list[str], similarity: float, n_threads: int = 1) -> pd.DataFrame:
    """Use the ChEMBL API to get similar compounds to a list of SMILES. Though multiple threads
    can be used, the rate limit is set to 5 calls per second not to overload the API.

    Args:
        smiles: list of SMILES strings to find similar molecules to.
        similarity: similarity threshold to use for the search. Value should be between 40 and 100. Defaults to 80.
        n_threads: Number of threads to use for searching the similar compounds. Defaults to 1.

    Returns:
        pd.DataFrame: a DataFrame with the similar molecules.
    """

    extracted = []

    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for smi in smiles:
            futures.append(executor.submit(get_similarity_compound_table, smi, similarity))
        for future in as_completed(futures):
            extracted.append(future.result())

    return pd.concat(extracted, ignore_index=True)
