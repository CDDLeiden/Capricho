"""Module for stereochemistry related functions"""

from typing import List, Union

from rdkit import Chem


def find_undefined_stereocenters(input: Union[Chem.Mol | str]) -> List[int]:
    """
    Find atoms that are stereocenters but have undefined chirality.

    Args:
        mol: Input RDKit molecule

    Returns:
        List[int]: List of atom indices that are undefined stereocenters
    """
    if input is None:
        return []

    if isinstance(input, str):
        mol = Chem.MolFromSmiles(input)

    elif isinstance(input, Chem.Mol):
        mol = input

    try:
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    except TypeError as e:
        raise TypeError(f"Something wront with input: {input}") from e

    undefined_stereo = []  # find atoms with undefined stereochemistry
    for atom_idx, chirality in chiral_centers:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            undefined_stereo.append(atom_idx)

    return undefined_stereo
