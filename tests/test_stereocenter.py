"""Test function to check for undefined stereocenters in a molecule."""

import unittest
from typing import List

from rdkit import Chem

from CompoundMapper.core.stereo import find_undefined_stereocenters


def has_undefined_stereocenter(mol: Chem.Mol) -> bool:
    """
    Check if molecule has any stereocenters with undefined chirality.

    Args:
        mol: Input RDKit molecule

    Returns:
        bool: True if molecule has undefined stereocenters, False otherwise
    """
    return len(find_undefined_stereocenters(mol)) > 0


class TestStereocenters(unittest.TestCase):
    def test_no_stereocenter(self):
        # Ethane has no stereocenters
        mol = Chem.MolFromSmiles("CC")
        self.assertFalse(has_undefined_stereocenter(mol))
        self.assertEqual(find_undefined_stereocenters(mol), [])

    def test_defined_stereocenter(self):
        # 2-butanol with defined stereochemistry
        mol = Chem.MolFromSmiles("CC[C@H](C)O")
        self.assertFalse(has_undefined_stereocenter(mol))
        self.assertEqual(find_undefined_stereocenters(mol), [])

    def test_undefined_stereocenter(self):
        # 2-butanol without defined stereochemistry
        mol = Chem.MolFromSmiles("CCC(C)O")
        self.assertTrue(has_undefined_stereocenter(mol))
        self.assertEqual(len(find_undefined_stereocenters(mol)), 1)

    def test_multiple_stereocenters(self):
        # 2,3-butanediol with one defined and one undefined stereocenter
        mol = Chem.MolFromSmiles("C[C@H](O)C(C)O")
        self.assertTrue(has_undefined_stereocenter(mol))
        undefined = find_undefined_stereocenters(mol)
        self.assertEqual(len(undefined), 1)

    def test_all_undefined_stereocenters(self):
        # 2,3-butanediol with both stereocenters undefined
        mol = Chem.MolFromSmiles("CC(O)C(C)O")
        self.assertTrue(has_undefined_stereocenter(mol))
        undefined = find_undefined_stereocenters(mol)
        self.assertEqual(len(undefined), 2)

    def test_none_input(self):
        self.assertFalse(has_undefined_stereocenter(None))
        self.assertEqual(find_undefined_stereocenters(None), [])

    def test_aromatic_compound(self):
        # Phenethyl alcohol has no stereocenters
        mol = Chem.MolFromSmiles("OCCc1ccccc1")
        self.assertFalse(has_undefined_stereocenter(mol))
        self.assertEqual(find_undefined_stereocenters(mol), [])

    def test_complex_molecule(self):
        # A more complex molecule with multiple stereocenters
        mol = Chem.MolFromSmiles("CC(C)[C@H](C(=O)N[C@@H](Cc1ccccc1)C(C)C)NC(=O)C(C)C")
        self.assertFalse(has_undefined_stereocenter(mol))
        self.assertEqual(find_undefined_stereocenters(mol), [])
