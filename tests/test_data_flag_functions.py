"""Tests for data_flag_functions module."""

import unittest

import pandas as pd

from Capricho.chembl.data_flag_functions import (
    flag_censored_activity_comment,
    flag_incompatible_units,
    flag_insufficient_assay_overlap,
    flag_inter_document_duplication,
    flag_missing_document_date,
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


class TestFlagMissingDocumentDate(unittest.TestCase):
    """Tests for flag_missing_document_date function."""

    def test_flag_missing_year(self):
        """Test that activities with missing year are flagged in processing comment."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
                "year": [2020, None, 2021],
                "pchembl_value": [6.0, 6.5, 7.0],
                "data_dropping_comment": [None, None, None],
            }
        )

        result = flag_missing_document_date(df)

        # Check that only row 1 (with None year) is flagged
        self.assertFalse(
            "Missing document date" in str(result.loc[0, "data_dropping_comment"]),
            "Row 0 should NOT be flagged (has year)",
        )
        self.assertTrue(
            "Missing document date" in str(result.loc[1, "data_dropping_comment"]),
            "Row 1 should be flagged (missing year)",
        )
        self.assertFalse(
            "Missing document date" in str(result.loc[2, "data_dropping_comment"]),
            "Row 2 should NOT be flagged (has year)",
        )

    def test_all_have_year(self):
        """Test that when all activities have year, nothing is flagged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
                "year": [2019, 2020, 2021],
                "pchembl_value": [6.0, 6.5, 7.0],
                "data_processing_comment": [None, None, None],
            }
        )

        result = flag_missing_document_date(df)

        # No rows should have "Missing document date" in processing comment
        for idx in range(len(result)):
            comment = result.loc[idx, "data_processing_comment"]
            self.assertFalse(
                comment and "Missing document date" in str(comment),
                f"Row {idx} should NOT be flagged (has year)",
            )

    def test_all_missing_year(self):
        """Test that when all activities lack year, all are flagged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
                "year": [None, None, None],
                "pchembl_value": [6.0, 6.5, 7.0],
                "data_dropping_comment": [None, None, None],
            }
        )

        result = flag_missing_document_date(df)

        # All rows should be flagged
        for idx in range(len(result)):
            self.assertTrue(
                "Missing document date" in str(result.loc[idx, "data_dropping_comment"]),
                f"Row {idx} should be flagged (missing year)",
            )

    def test_missing_year_column(self):
        """Test that function handles missing year column gracefully."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
                "pchembl_value": [6.0, 6.5],
                "data_processing_comment": [None, None],
            }
        )

        result = flag_missing_document_date(df)

        # Should return DataFrame unchanged
        self.assertEqual(len(result), 2)
        self.assertTrue("year" not in result.columns)

    def test_preserves_existing_comments(self):
        """Test that existing processing comments are preserved when flagging."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
                "year": [None, None],
                "pchembl_value": [6.0, 6.5],
                "data_dropping_comment": ["Existing comment", None],
            }
        )

        result = flag_missing_document_date(df)

        # Check that existing comment is preserved and new flag is added
        self.assertIn("Existing comment", result.loc[0, "data_dropping_comment"])
        self.assertIn("Missing document date", result.loc[0, "data_dropping_comment"])
        self.assertIn("Missing document date", result.loc[1, "data_dropping_comment"])


class TestFlagIncompatibleUnits(unittest.TestCase):
    """Tests for flag_incompatible_units function."""

    def test_flag_incompatible_units(self):
        """Test that activities with non-convertible units are flagged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3", "CHEMBL4"],
                "standard_units": ["nM", "%", "µM", "ug.mL-1"],
                "standard_value": [100.0, 50.0, 10.0, 5.0],
                "pchembl_value": [7.0, None, 8.0, None],
                "data_dropping_comment": [None, None, None, None],
            }
        )

        result = flag_incompatible_units(df)

        # Check that rows with "%" and "ug.mL-1" are flagged, but "nM" and "µM" are not
        self.assertFalse(
            "Incompatible units" in str(result.loc[0, "data_dropping_comment"]),
            "Row 0 (nM) should NOT be flagged",
        )
        self.assertTrue(
            "Incompatible units" in str(result.loc[1, "data_dropping_comment"]),
            "Row 1 (%) should be flagged",
        )
        self.assertFalse(
            "Incompatible units" in str(result.loc[2, "data_dropping_comment"]),
            "Row 2 (µM) should NOT be flagged",
        )
        self.assertTrue(
            "Incompatible units" in str(result.loc[3, "data_dropping_comment"]),
            "Row 3 (L) should be flagged",
        )

    def test_all_compatible_units(self):
        """Test that when all units are compatible, nothing is flagged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2", "CHEMBL3"],
                "standard_units": ["nM", "uM", "mM"],
                "standard_value": [100.0, 10.0, 1.0],
                "pchembl_value": [7.0, 8.0, 6.0],
                "data_dropping_comment": [None, None, None],
            }
        )

        result = flag_incompatible_units(df)

        # No rows should be flagged
        for idx in range(len(result)):
            comment = result.loc[idx, "data_dropping_comment"]
            self.assertFalse(
                comment and "Incompatible units" in str(comment),
                f"Row {idx} should NOT be flagged (compatible unit)",
            )

    def test_null_units_not_flagged(self):
        """Test that null/NA units are not flagged (they're handled separately)."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL2"],
                "standard_units": [None, "nM"],
                "standard_value": [7.5, 100.0],
                "pchembl_value": [7.5, 7.0],
                "data_dropping_comment": [None, None],
            }
        )

        result = flag_incompatible_units(df)

        # Null units should not be flagged
        self.assertFalse(
            result.loc[0, "data_dropping_comment"]
            and "Incompatible units" in str(result.loc[0, "data_dropping_comment"]),
            "Row 0 (null unit) should NOT be flagged",
        )

    def test_missing_standard_units_column(self):
        """Test that function handles missing standard_units column gracefully."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["CHEMBL1"],
                "pchembl_value": [6.0],
                "data_dropping_comment": [None],
            }
        )

        result = flag_incompatible_units(df)

        # Should return DataFrame unchanged
        self.assertEqual(len(result), 1)


class TestFlagInsufficientAssayOverlap(unittest.TestCase):
    """Tests for flag_insufficient_assay_overlap function with and without metadata matching."""

    def test_basic_overlap_filtering(self):
        """Test that assays without sufficient overlap are flagged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["MOL1", "MOL2", "MOL3", "MOL4", "MOL5"],
                "assay_chembl_id": ["ASSAY1", "ASSAY1", "ASSAY2", "ASSAY2", "ASSAY2"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1", "TARGET1", "TARGET1"],
                "document_chembl_id": ["DOC1", "DOC1", "DOC1", "DOC1", "DOC1"],
                "pchembl_value": [6.0, 6.5, 7.0, 7.5, 8.0],
                "data_dropping_comment": [None, None, None, None, None],
            }
        )

        # ASSAY1 has MOL1, MOL2 (2 compounds)
        # ASSAY2 has MOL3, MOL4, MOL5 (3 compounds)
        # Overlap between them: 0 compounds
        # With min_overlap=2, both assays should be flagged

        result = flag_insufficient_assay_overlap(df, min_overlap=2)

        # All activities should be flagged since no assay pair meets min_overlap
        flagged_mask = result["data_dropping_comment"].str.contains("Insufficient assay overlap", na=False)
        self.assertTrue(flagged_mask.all(), "All activities should be flagged when overlap is insufficient")

    def test_sufficient_overlap_not_flagged(self):
        """Test that assays with sufficient overlap are not flagged."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["MOL1", "MOL2", "MOL1", "MOL2", "MOL3"],
                "assay_chembl_id": ["ASSAY1", "ASSAY1", "ASSAY2", "ASSAY2", "ASSAY2"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value": [6.0, 6.5, 6.1, 6.6, 7.0],
                "data_dropping_comment": [None, None, None, None, None],
            }
        )

        # ASSAY1 has MOL1, MOL2
        # ASSAY2 has MOL1, MOL2, MOL3
        # Overlap: 2 compounds (MOL1, MOL2)
        # With min_overlap=2, both assays should NOT be flagged

        result = flag_insufficient_assay_overlap(df, min_overlap=2)

        # No activities should be flagged
        flagged_mask = result["data_dropping_comment"].str.contains("Insufficient assay overlap", na=False)
        self.assertFalse(flagged_mask.any(), "No activities should be flagged when overlap is sufficient")

    def test_multiple_targets_handled_independently(self):
        """Test that overlap checking is done independently for each target."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["MOL1", "MOL2", "MOL1", "MOL2"],
                "assay_chembl_id": ["ASSAY1", "ASSAY1", "ASSAY2", "ASSAY2"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET2", "TARGET2"],  # Different targets
                "pchembl_value": [6.0, 6.5, 7.0, 7.5],
                "data_dropping_comment": [None, None, None, None],
            }
        )

        # TARGET1: ASSAY1 has MOL1, MOL2 (only one assay, can't form pair)
        # TARGET2: ASSAY2 has MOL1, MOL2 (only one assay, can't form pair)
        # With min_overlap=2, nothing should be flagged (single assays per target)

        result = flag_insufficient_assay_overlap(df, min_overlap=2)

        # No activities should be flagged (no pairs to compare)
        flagged_mask = result["data_dropping_comment"].str.contains("Insufficient assay overlap", na=False)
        self.assertFalse(
            flagged_mask.any(), "No activities should be flagged when each target has only one assay"
        )

    def test_min_overlap_zero_skips_filtering(self):
        """Test that min_overlap=0 skips the filtering entirely."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["MOL1", "MOL2"],
                "assay_chembl_id": ["ASSAY1", "ASSAY2"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "pchembl_value": [6.0, 6.5],
                "data_dropping_comment": [None, None],
            }
        )

        result = flag_insufficient_assay_overlap(df, min_overlap=0)

        # No activities should be flagged
        flagged_mask = result["data_dropping_comment"].str.contains("Insufficient assay overlap", na=False)
        self.assertFalse(flagged_mask.any(), "min_overlap=0 should skip all filtering")

    def test_skips_size_flagged_assays_goldilocks(self):
        """Test that assays already flagged for size issues are skipped in overlap checking (goldilocks approach)."""
        df = pd.DataFrame(
            {
                "molecule_chembl_id": ["MOL1", "MOL2", "MOL3", "MOL4", "MOL5", "MOL6"],
                "assay_chembl_id": ["ASSAY1", "ASSAY1", "ASSAY2", "ASSAY2", "ASSAY3", "ASSAY3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1", "TARGET1", "TARGET1", "TARGET1"],
                "document_chembl_id": ["DOC1", "DOC1", "DOC2", "DOC2", "DOC3", "DOC3"],
                "pchembl_value": [6.0, 6.5, 6.0, 7.0, 6.0, 8.0],
                "data_dropping_comment": [
                    "Assay size < 20",
                    "Assay size < 20",
                    None,
                    None,
                    None,
                    None,
                ],
            }
        )

        # ASSAY1 is already flagged for size (has MOL1, MOL2)
        # ASSAY2 has MOL3, MOL4
        # ASSAY3 has MOL5, MOL6
        # ASSAY1 shares MOL1 with no one, but it's already size-flagged so should be skipped
        # ASSAY2 and ASSAY3 have no overlap, so both should be flagged for insufficient overlap

        result = flag_insufficient_assay_overlap(df, min_overlap=1)

        # ASSAY1 should still only have size flag, not overlap flag
        assay1_comments = result[result["assay_chembl_id"] == "ASSAY1"]["data_dropping_comment"]
        for comment in assay1_comments:
            self.assertIn("Assay size < 20", comment)
            self.assertNotIn("Insufficient assay overlap", comment)

        # ASSAY2 and ASSAY3 should have overlap flag (they don't overlap with each other or ASSAY1)
        assay2_comments = result[result["assay_chembl_id"] == "ASSAY2"]["data_dropping_comment"]
        for comment in assay2_comments:
            self.assertIn("Insufficient assay overlap", comment)

        assay3_comments = result[result["assay_chembl_id"] == "ASSAY3"]["data_dropping_comment"]
        for comment in assay3_comments:
            self.assertIn("Insufficient assay overlap", comment)


if __name__ == "__main__":
    unittest.main()
