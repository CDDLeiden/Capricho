"""Utility functions for further processing the standardized smiles by removing mixtures and salts."""

import numpy as np
import re

MIXTURE_REGEX = re.compile(
    r"^("  # start of the string
    r"Na\+?|Cl\-?|Br\-?|K\+?|F\-?|I\-?|Ca\+{2}?|Mg\+{2}?|Zn\+{2}?|OH\-?|"  # Salts
    r"\[Na\+\]?|\[K\+\]?|\[Cl-\]?|\[Br-\]?|\[I-\]?|\[Zn\+2\]?|"  # Salts
    r"CCCC\(=O\)\[O-\]|"  # Butyrate
    r"CCCCC\(=O\)\[O-\]|"  # Pentanoate
    r"C1=CC=C\(C=C1\)C\(=O\)\[O-\]|"  # Benzoate
    r"\[O-\]\[Cl\+3\]\(\[O-\]\)\(\[O-\]\)\[O-\]|"  # Perchlorate
    r"C1=CC=NC=C1|"  # Pyridine
    r"O=CN\(C\)C|"  # N,N-Dimethylformamide
    r"\[NO\]|"  # Nitric oxide, corrected
    r"c1ccc\(\[B-\]\(c2ccccc2\)\(c2ccccc2\)c2ccccc2\)cc1|"  # Tetraphenylborate
    r"^O$|^N$"  # Match full strings for O and N
    r")$"  # end of the string
)


def clean_mixtures(smi):
    """Removes mixtures/salts from a SMILES tring that aren't captured by the
    `chembl_structure_pipeline` using MIXTURE_REGEX.

    Args:
        smi: smiles string

    Returns:
        str: smiles string without mixture.
    """
    smiles = np.unique(smi.split(".")).tolist()
    if len(smiles) == 1:
        return smiles[0]
    else:
        no_salt_smiles = [MIXTURE_REGEX.sub("", smi) for smi in smiles]
        no_salt_smiles = [smi for smi in no_salt_smiles if smi]
        if len(no_salt_smiles) == 0:
            return "."  # in this case wil drop the smiles
        return ".".join(no_salt_smiles)
