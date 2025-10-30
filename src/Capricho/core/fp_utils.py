"""Utility functions to calculate fingerprints & identify molecules to be treated as identical"""

from functools import partial

import numpy as np
from job_tqdflex import ParallelApplier
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


def smi_to_morganFP(smi, radius: int = 2, nBits=2048, useChirality=False, **kwargs) -> np.ndarray:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Invalid SMILES detected: {smi}")
        return None
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=nBits, includeChirality=useChirality, **kwargs
    )
    return morgan_gen.GetFingerprintAsNumPy(mol).reshape(1, -1)


def smi_to_RDKitFP(smi, minPath=1, maxPath=7, nBits=2048, **kwargs) -> np.ndarray:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"Invalid SMILES detected: {smi}")
        return None
    rdkit_gen = rdFingerprintGenerator.GetRDKitFPGenerator(
        minPath=minPath, maxPath=maxPath, fpSize=nBits, **kwargs
    )
    return rdkit_gen.GetFingerprintAsNumPy(mol).reshape(1, -1)


def calculate_mixed_FPs(
    smiles: list,
    n_jobs: int = 8,
    morgan_kwargs: dict = None,
    rdkit_kwargs: dict = None,
    return_stacked: bool = False,
    chunk_size: int = 50,
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
        chunk_size (int): chunk size to use for the parallel applier. Defaults to 50.

    Returns:
        np.ndarray: a mixed fingerprint for the input smiles
    """
    if morgan_kwargs is None:
        morgan_kwargs = {}
    if rdkit_kwargs is None:
        rdkit_kwargs = {}

    morganfunc = partial(smi_to_morganFP, **morgan_kwargs)
    rdkitfpfunc = partial(smi_to_RDKitFP, **rdkit_kwargs)

    morgan_applier = ParallelApplier(
        morganfunc,
        smiles,
        n_jobs=n_jobs,
        backend="loky",
        show_progress=True,
        chunk_size=chunk_size,
        custom_desc="Calculating Morgan FPs",
    )
    rdkit_applier = ParallelApplier(
        rdkitfpfunc,
        smiles,
        n_jobs=n_jobs,
        backend="loky",
        show_progress=True,
        chunk_size=chunk_size,
        custom_desc="Calculating RDKit FPs",
    )

    morgan_fps = morgan_applier()
    rdkit_fps = rdkit_applier()

    if return_stacked:
        return np.concatenate([np.concatenate(morgan_fps), np.concatenate(rdkit_fps)], axis=1).shape
    else:
        return [np.concatenate([morfp, rdkfp], axis=1) for morfp, rdkfp in zip(morgan_fps, rdkit_fps)]
