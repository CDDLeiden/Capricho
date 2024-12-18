import numpy as np
import pandas as pd
from chemFilters.chem.standardizers import ChemStandardizer

from ..core.fp_utils import calculate_mixed_FPs
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import repeated_indices_from_array_series
from .core import get_multiple_compounds


def get_and_curate_multiple_compounds_result(
    cpd_list: list[str],
    input_type: str = "name",
    chirality: bool = True,
    n_jobs: int = 1,
    n_threads: int = 4,
) -> pd.DataFrame:
    """
    Will fetch a list of inputs from PubChem, treating them as the input type specified.
    For the list of compounds fetched, it will standardize the SMILES, remove mixtures,
    and calculate fingerprints to identify repeated compounds.

    Same compounds will have the same compound identifier in the "repeats" column.

    Args:
        cpd_list: List of compound identifiers to fetch from PubChem
        input_type: Type of input. Available: name, smiles, sdf, inchi, inchikey, formula
        chirality: Use chirality to identify if compounds are the same with fingerprints
        n_jobs: Number of parallel jobs for calculating fingerprints. 1 is advised. Defaults to 1.
        n_threads: Number of parallel threads for fetching compounds. 4 is advised. Defaults to 4.

    Returns:
        df: DataFrame with curated results
    """
    pubchempy_results = get_multiple_compounds(cpd_list, input_type, n_jobs=n_threads)

    stdzer = ChemStandardizer(from_smi=True, n_jobs=n_jobs, verbose=False)
    curated = []

    for inp, res in zip(cpd_list, pubchempy_results):
        empty_entry = pd.DataFrame(
            [
                {
                    input_type: inp,
                    "smiles": np.nan,
                    "pubchem_cid": np.nan,
                    "inchi": np.nan,
                    "inchikey": np.nan,
                    "isomeric_smiles": np.nan,
                    "canonical_smiles": np.nan,
                    "iupac_name": np.nan,
                    "synonyms": np.nan,
                }
            ]
        )

        if not res:
            curated.append(empty_entry)
            continue

        to_curate = [
            {
                input_type: inp,
                "smiles": (r.isomeric_smiles if chirality else r.canonical_smiles),
                "pubchem_cid": r.cid,
                "inchi": r.inchi,
                "inchikey": r.inchikey,
                "isomeric_smiles": r.isomeric_smiles,
                "canonical_smiles": r.canonical_smiles,
                "iupac_name": r.iupac_name,
                "synonyms": r.synonyms,
            }
            for r in res
        ]

        df = (
            pd.DataFrame(to_curate)
            .assign(standard_smiles=lambda x: stdzer(x["smiles"]))
            .dropna(subset=["standard_smiles"])  # drop if no structure is found
            .query("standard_smiles.notna()")
            .assign(final_smiles=lambda x: x["standard_smiles"].apply(clean_mixtures))
            .drop(columns=["standard_smiles", "smiles"])
            .rename(columns={"final_smiles": "smiles"})
        )

        if df.shape[0] == 1:
            curated.append(df)
            continue
        elif df.shape[0] == 0:
            curated.append(empty_entry)
            continue
        else:
            fps = calculate_mixed_FPs(
                df["smiles"].tolist(), n_jobs=n_jobs, morgan_kwargs={"useChirality": chirality}
            )
            df = df.assign(fps=fps).assign(repeats=False)
            repeats_idxs = repeated_indices_from_array_series(df["fps"])
            for repeats in repeats_idxs:
                df.loc[repeats, "repeats"] = df.loc[repeats, input_type]

        curated.append(df)

    curated = (
        pd.concat(curated, ignore_index=True)
        .drop(columns=["fps"])
        .assign(repeats=lambda x: x["repeats"].fillna(False))
    )
    return curated
