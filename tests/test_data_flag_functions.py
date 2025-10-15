"""Tests for data_flag_functions module."""

import unittest

import pandas as pd

from Capricho.chembl.data_flag_functions import (
    flag_censored_activity_comment,
    flag_inter_document_duplication,
)


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


class TestFlagInterDocumentDuplication(unittest.TestCase):
    """Validates censored activity comment handling and inter-document duplication detection for discrete measurements."""

    def test_only_flags_discrete_measurements(self):
        """Test that only discrete measurements (standard_relation='=') are flagged as duplicates."""
        test_data = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"] * 6,
                "standard_smiles": ["CCO"] * 6,
                "canonical_smiles": ["CCO"] * 6,
                "pchembl_value": [6.0, 6.0, 6.0, 6.0, 5.0, 5.0],
                "standard_relation": ["=", "=", "<", "<", "=", "="],
                "target_chembl_id": ["TARGET1"] * 6,
                "mutation": ["WT"] * 6,
                "target_organism": ["Homo sapiens"] * 6,
                "document_chembl_id": ["DOC1", "DOC2", "DOC3", "DOC4", "DOC5", "DOC6"],
                "data_processing_comment": [""] * 6,
            }
        )

        result = flag_inter_document_duplication(test_data)

        # Only rows 0 and 1 (discrete measurements with pchembl=6.0, relation='=') should be flagged
        # Rows 2 and 3 (censored measurements with pchembl=6.0, relation='<') should NOT be flagged
        # Rows 4 and 5 (discrete measurements with pchembl=5.0, relation='=') should be flagged
        flagged_mask = result["data_processing_comment"].str.contains(
            "pChEMBL Duplication Across Documents", na=False
        )

        # Check that rows 0, 1, 4, 5 are flagged (discrete measurements)
        self.assertTrue(flagged_mask.iloc[0], "Row 0 should be flagged (discrete, pchembl=6.0)")
        self.assertTrue(flagged_mask.iloc[1], "Row 1 should be flagged (discrete, pchembl=6.0)")
        self.assertTrue(flagged_mask.iloc[4], "Row 4 should be flagged (discrete, pchembl=5.0)")
        self.assertTrue(flagged_mask.iloc[5], "Row 5 should be flagged (discrete, pchembl=5.0)")

        # Check that rows 2 and 3 are NOT flagged (censored measurements)
        self.assertFalse(flagged_mask.iloc[2], "Row 2 should NOT be flagged (censored, relation='<')")
        self.assertFalse(flagged_mask.iloc[3], "Row 3 should NOT be flagged (censored, relation='<')")

    def test_no_discrete_measurements(self):
        """Test that when all measurements are censored, nothing is flagged."""
        test_data = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"] * 4,
                "standard_smiles": ["CCO"] * 4,
                "canonical_smiles": ["CCO"] * 4,
                "pchembl_value": [6.0, 6.0, 5.0, 5.0],
                "standard_relation": ["<", "<", ">", ">"],
                "target_chembl_id": ["TARGET1"] * 4,
                "mutation": ["WT"] * 4,
                "target_organism": ["Homo sapiens"] * 4,
                "document_chembl_id": ["DOC1", "DOC2", "DOC3", "DOC4"],
                "data_processing_comment": [""] * 4,
            }
        )

        result = flag_inter_document_duplication(test_data)

        # No rows should be flagged since all are censored
        flagged_mask = result["data_processing_comment"].str.contains(
            "pChEMBL Duplication Across Documents", na=False
        )
        self.assertFalse(flagged_mask.any(), "No censored measurements should be flagged")

    def test_missing_standard_relation_column(self):
        """Test that function handles missing standard_relation column gracefully."""
        test_data = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"] * 2,
                "standard_smiles": ["CCO"] * 2,
                "canonical_smiles": ["CCO"] * 2,
                "pchembl_value": [6.0, 6.0],
                "target_chembl_id": ["TARGET1"] * 2,
                "mutation": ["WT"] * 2,
                "target_organism": ["Homo sapiens"] * 2,
                "document_chembl_id": ["DOC1", "DOC2"],
                "data_processing_comment": [""] * 2,
            }
        )

        result = flag_inter_document_duplication(test_data)

        # Function should return dataframe unchanged
        self.assertEqual(len(result), 2)
        self.assertTrue(result["data_processing_comment"].str.strip().eq("").all())


if __name__ == "__main__":
    unittest.main()
