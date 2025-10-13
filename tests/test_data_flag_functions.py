"""ABOUTME: Tests for data_flag_functions module.
ABOUTME: Validates that censored activity comments are correctly identified and standard_relation is changed.
"""

import unittest

import pandas as pd

from Capricho.chembl.data_flag_functions import flag_censored_activity_comment


class TestFlagCensoredActivityComment(unittest.TestCase):
    def test_flag_inconclusive_comment(self):
        """Test that 'Inconclusive' activity_comment with '=' relation is corrected to '<'."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [6.0],
                "standard_relation": ["="],
                "activity_comment": ["Inconclusive"],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        self.assertEqual(result.loc[0, "standard_relation"], "<")
        self.assertIn("Corrected standard_relation", str(result.loc[0, "data_processing_comment"]))

    def test_flag_not_active_comment(self):
        """Test that 'Not Active' activity_comment with '=' relation is corrected to '<'."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [5.5],
                "standard_relation": ["="],
                "activity_comment": ["Not Active"],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        self.assertEqual(result.loc[0, "standard_relation"], "<")
        self.assertIn("Corrected standard_relation", str(result.loc[0, "data_processing_comment"]))

    def test_flag_inactive_comment(self):
        """Test that 'inactive at 10 uM' activity_comment with '=' relation is corrected to '<'."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [5.0],
                "standard_relation": ["="],
                "activity_comment": ["inactive at 10 uM"],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        self.assertEqual(result.loc[0, "standard_relation"], "<")
        self.assertIn("Corrected standard_relation", str(result.loc[0, "data_processing_comment"]))

    def test_no_change_for_active_comment(self):
        """Test that 'Active' activity_comment with '=' relation remains unchanged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [7.0],
                "standard_relation": ["="],
                "activity_comment": ["Active"],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        self.assertEqual(result.loc[0, "standard_relation"], "=")
        self.assertTrue(
            pd.isna(result.loc[0, "data_processing_comment"])
            or result.loc[0, "data_processing_comment"] in [None, "", pd.NA]
        )

    def test_no_change_for_already_correct_relation(self):
        """Test that 'Inactive' activity_comment with '>' relation remains unchanged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [6.5],
                "standard_relation": [">"],
                "activity_comment": ["Inactive"],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        self.assertEqual(result.loc[0, "standard_relation"], ">")

    def test_no_change_for_null_comment(self):
        """Test that null activity_comment with '=' relation remains unchanged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [6.0],
                "standard_relation": ["="],
                "activity_comment": [None],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        self.assertEqual(result.loc[0, "standard_relation"], "=")

    def test_case_insensitive_matching(self):
        """Test that keyword matching is case-insensitive."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
                "pchembl_value": [6.0, 5.5, 5.0],
                "standard_relation": ["=", "=", "="],
                "activity_comment": ["INCONCLUSIVE", "not active", "InAcTiVe"],
                "data_processing_comment": [None, None, None],
            }
        )

        result = flag_censored_activity_comment(df)

        self.assertEqual(result.loc[0, "standard_relation"], "<")
        self.assertEqual(result.loc[1, "standard_relation"], "<")
        self.assertEqual(result.loc[2, "standard_relation"], "<")

    def test_batch_correction(self):
        """Test that multiple rows are corrected in a single call."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3", "CHEMBL4", "CHEMBL5", "CHEMBL6"],
                "pchembl_value": [6.0, 5.5, 7.0, 6.5, 5.0, 6.0],
                "standard_relation": ["=", "=", "=", ">", "=", "="],
                "activity_comment": [
                    "Inconclusive",
                    "Not Active",
                    "Active",
                    "Inactive",
                    "inactive at 10 uM",
                    None,
                ],
                "data_processing_comment": [None, None, None, None, None, None],
            }
        )

        result = flag_censored_activity_comment(df)

        # Check each row
        self.assertEqual(
            result.loc[0, "standard_relation"], "<", "CHEMBL1: Should change '=' to '<' (Inconclusive)"
        )
        self.assertEqual(
            result.loc[1, "standard_relation"], "<", "CHEMBL2: Should change '=' to '<' (Not Active)"
        )
        self.assertEqual(result.loc[2, "standard_relation"], "=", "CHEMBL3: Should remain '=' (Active)")
        self.assertEqual(
            result.loc[3, "standard_relation"], ">", "CHEMBL4: Should remain '>' (already correct)"
        )
        self.assertEqual(
            result.loc[4, "standard_relation"], "<", "CHEMBL5: Should change '=' to '<' (inactive at 10 uM)"
        )
        self.assertEqual(
            result.loc[5, "standard_relation"], "=", "CHEMBL6: Should remain '=' (no activity_comment)"
        )

    def test_missing_activity_comment_column(self):
        """Test that function handles missing activity_comment column gracefully."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [6.0],
                "standard_relation": ["="],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        # Should return DataFrame unchanged
        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "standard_relation"], "=")

    def test_missing_standard_relation_column(self):
        """Test that function handles missing standard_relation column gracefully."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [6.0],
                "activity_comment": ["Inactive"],
                "data_processing_comment": [None],
            }
        )

        result = flag_censored_activity_comment(df)

        # Should return DataFrame unchanged
        self.assertEqual(len(result), 1)


if __name__ == "__main__":
    unittest.main()
