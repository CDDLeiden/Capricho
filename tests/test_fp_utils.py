import unittest

import numpy as np

from CompoundMapper.core import fp_utils


class TestFPUtils(unittest.TestCase):
    def setUp(self):
        self.test_smiles = ["CC", "CCO", "c1ccccc1"]

    def test_smi_to_morganFP(self):
        fp = fp_utils.smi_to_morganFP(self.test_smiles[0])
        self.assertIsInstance(fp, np.ndarray)
        self.assertEqual(fp.shape, (1, 2048))

    def test_smi_to_RDKitFP(self):
        fp = fp_utils.smi_to_RDKitFP(self.test_smiles[0])
        self.assertIsInstance(fp, np.ndarray)
        self.assertEqual(fp.shape, (1, 2048))

    def test_calculate_mixed_FPs(self):
        fps = fp_utils.calculate_mixed_FPs(self.test_smiles, n_jobs=1)
        self.assertEqual(len(fps), len(self.test_smiles))
        self.assertIsInstance(fps[0], np.ndarray)
        self.assertEqual(fps[0].shape, (1, 4096))  # 2048 (Morgan) + 2048 (RDKit)

    def test_calculate_mixed_FPs_parallel(self):
        fps = fp_utils.calculate_mixed_FPs(self.test_smiles, n_jobs=2)
        self.assertEqual(len(fps), len(self.test_smiles))
        self.assertIsInstance(fps[0], np.ndarray)
        self.assertEqual(fps[0].shape, (1, 4096))


if __name__ == "__main__":
    unittest.main()
