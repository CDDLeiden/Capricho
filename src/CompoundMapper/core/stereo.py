"""Module for stereochemistry related functions"""

from typing import List

from rdkit import Chem


def find_undefined_stereocenters(mol: Chem.Mol) -> List[int]:
    """
    Find atoms that are stereocenters but have undefined chirality.

    Args:
        mol: Input RDKit molecule

    Returns:
        List[int]: List of atom indices that are undefined stereocenters
    """
    if mol is None:
        return []

    # Get all potential stereocenters
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

    # Find atoms with unassigned stereochemistry
    undefined_stereo = []
    for atom_idx, chirality in chiral_centers:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            undefined_stereo.append(atom_idx)

    return undefined_stereo
