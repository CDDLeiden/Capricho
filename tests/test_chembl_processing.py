import unittest

import pandas as pd

from CompoundMapper.chembl import processing


class TestChemblProcessing(unittest.TestCase):
    def setUp(self):
        self.sample_df = pd.DataFrame(
            {
                "standard_value": [1, 0.2, 1],
                "standard_units": ["nM", "µM", "uM"],
                "expected_pchembl_value": [9.0, 6.7, 6.5],
                "pchembl_value": [None, None, None],
                "data_validity_comment": [None, None, "Outside typical range"],
                "potential_duplicate": [0, 0, 1],
                "assay_type": ["B", "F", "B"],
                "standard_type": ["IC50", "Ki", "EC50"],
            }
        )

    def test_convert_to_log10(self):
        result = processing.convert_to_log10(self.sample_df)
        self.assertAlmostEqual(result["pchembl_value"].iloc[0], 9.0, places=2)
        self.assertAlmostEqual(result["pchembl_value"].iloc[1], 6.7, places=2)
        self.assertAlmostEqual(result["pchembl_value"].iloc[2], 6, places=2)
        # Raise an error when pChEMBL value is already present...
        self.assertRaises(ValueError, processing.convert_to_log10, self.sample_df.assign(pchembl_value=1))
        self.assertRaises(
            ValueError, processing.convert_to_log10, self.sample_df.drop(columns=["pchembl_value"])
        )

    def test_process_bioactivities(self):
        result = processing.process_bioactivities(self.sample_df)
        self.assertEqual(len(result), 2)  # One row should be filtered out
        self.assertNotIn("data_validity_description", result.columns)
        self.assertNotIn("potential_duplicate", result.columns)


if __name__ == "__main__":
    unittest.main()
