"""Utility functions to calculate fingerprints & identify molecules to be treated as identical"""

from functools import partial

import numpy as np
from joblib import Parallel, delayed, parallel_config
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdmolops
from tqdm import tqdm


def smi_to_morganFP(smi, radius: int = 2, nBits=2048, useChirality=False, **kwargs) -> np.ndarray:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Invalid SMILES detected: {smi}")
        return None
    numpy_fp = np.zeros((1, nBits), dtype=bool)
    morgan_fp = AllChem.GetMorganFingerprintAsBitVect(
        mol, radius, nBits=nBits, useChirality=useChirality, **kwargs
    )
    DataStructs.ConvertToNumpyArray(morgan_fp, numpy_fp)
    return numpy_fp.reshape(1, -1)


def smi_to_RDKitFP(smi, minPath=1, maxPath=7, nBits=2048, **kwargs) -> np.ndarray:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Invalid SMILES detected: {smi}")
        return None
    numpy_fp = np.zeros((1, nBits), dtype=bool)
    morgan_fp = rdmolops.RDKFingerprint(mol, minPath=minPath, maxPath=maxPath, fpSize=nBits, **kwargs)
    DataStructs.ConvertToNumpyArray(morgan_fp, numpy_fp)
    return numpy_fp.reshape(1, -1)


def calculate_mixed_FPs(
    smiles: list,
    n_jobs: int = 8,
    morgan_kwargs: dict = None,
    rdkit_kwargs: dict = None,
    return_stacked: bool = False,
):
    """Outputs a mixed fingerprint used for compound identification. The motivation for this is
    that either of the fingerprints can fail to identify the same compound, but the combination
    is less prone to failure in this regard.

    Args:
        smiles (list): a list of smiles for which to compoute the mixed fingerprint
        n_jobs (int): number of jobs to run the fp calculation in parallel. Defaults to 1.
        morgan_kwargs (dict): keyword arguments for the morgan fingerprints. Defaults to None.
        rdkit_kwargs (dict): keyword arguments for the rdkit path fingerprints. Defaults to None.
        return_stacked (bool): if true, will return the stacked fingerprints instead of a list of
            numpy arrays. Defaults to False.

    Returns:
        np.ndarray: a mixed fingerprint for the input smiles
    """
    if morgan_kwargs is None:
        morgan_kwargs = {}
    if rdkit_kwargs is None:
        rdkit_kwargs = {}

    morganfunc = partial(smi_to_morganFP, **morgan_kwargs)
    rdkitfpfunc = partial(smi_to_RDKitFP, **rdkit_kwargs)

    with parallel_config(backend="loky", n_jobs=n_jobs):
        morgan_fps = Parallel()(
            delayed(morganfunc)(smi) for smi in tqdm(smiles, desc="Calculating Morgan fingerprints")
        )
        rdkit_fps = Parallel()(
            delayed(rdkitfpfunc)(smi) for smi in tqdm(smiles, desc="Calculating RDKit fingerprints")
        )
    if return_stacked:
        return np.concatenate([np.concatenate(morgan_fps), np.concatenate(rdkit_fps)], axis=1).shape
    else:
        return [np.concatenate([morfp, rdkfp], axis=1) for morfp, rdkfp in zip(morgan_fps, rdkit_fps)]
