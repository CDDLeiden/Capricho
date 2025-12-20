"""Tests for analysis module."""

import unittest

import pandas as pd

from Capricho.analysis import (
    DroppingComment,
    ProcessingComment,
    deaggregate_data,
    explode_assay_comparability,
    extract_assay_threshold,
    get_all_comments,
    normalize_comment_pattern,
)


class TestNormalizeCommentPattern(unittest.TestCase):
    """Tests for normalize_comment_pattern function."""

    def test_normalize_assay_size_min(self):
        """Test that assay size minimum patterns are normalized correctly."""
        self.assertEqual(normalize_comment_pattern("Assay size < 20"), "Assay size <")
        self.assertEqual(normalize_comment_pattern("Assay size < 10"), "Assay size <")
        self.assertEqual(normalize_comment_pattern("Assay size < 50"), "Assay size <")

    def test_normalize_assay_size_max(self):
        """Test that assay size maximum patterns are normalized correctly."""
        self.assertEqual(normalize_comment_pattern("Assay size > 100"), "Assay size >")
        self.assertEqual(normalize_comment_pattern("Assay size > 200"), "Assay size >")
        self.assertEqual(normalize_comment_pattern("Assay size > 50"), "Assay size >")

    def test_normalize_static_comments(self):
        """Test that static comments remain unchanged."""
        static_comments = [
            "Unit Annotation Error",
            "Potential Duplicate",
            "Data Validity Comment Present",
            "Undefined Stereochemistry",
            "Salt/solvent removed",
            "Calculated pChEMBL",
        ]
        for comment in static_comments:
            self.assertEqual(
                normalize_comment_pattern(comment), comment, f"Comment '{comment}' should not be normalized"
            )


class TestExtractAssayThreshold(unittest.TestCase):
    """Tests for extract_assay_threshold function."""

    def test_extract_min_threshold(self):
        """Test extraction of minimum assay size threshold."""
        self.assertEqual(extract_assay_threshold("Assay size < 20"), "20")
        self.assertEqual(extract_assay_threshold("Assay size < 10"), "10")
        self.assertEqual(extract_assay_threshold("Assay size < 50"), "50")

    def test_extract_max_threshold(self):
        """Test extraction of maximum assay size threshold."""
        self.assertEqual(extract_assay_threshold("Assay size > 100"), "100")
        self.assertEqual(extract_assay_threshold("Assay size > 200"), "200")
        self.assertEqual(extract_assay_threshold("Assay size > 500"), "500")

    def test_extract_non_assay_size_comment(self):
        """Test that non-assay-size comments return empty string."""
        self.assertEqual(extract_assay_threshold("Unit Annotation Error"), "")
        self.assertEqual(extract_assay_threshold("Potential Duplicate"), "")
        self.assertEqual(extract_assay_threshold(""), "")


class TestGetAllComments(unittest.TestCase):
    """Tests for get_all_comments function."""

    def test_returns_list(self):
        """Test that function returns a list."""
        result = get_all_comments()
        self.assertIsInstance(result, list)

    def test_contains_all_enum_values(self):
        """Test that all enum values are present in the returned list."""
        result = get_all_comments()

        # Check ProcessingComment values
        for comment in ProcessingComment:
            self.assertIn(
                comment.value, result, f"ProcessingComment '{comment.value}' should be in get_all_comments()"
            )

        # Check DroppingComment values
        for comment in DroppingComment:
            self.assertIn(
                comment.value, result, f"DroppingComment '{comment.value}' should be in get_all_comments()"
            )

    def test_no_duplicates(self):
        """Test that there are no duplicate values in the returned list."""
        result = get_all_comments()
        self.assertEqual(len(result), len(set(result)), "get_all_comments() should not contain duplicates")

    def test_assay_size_patterns_present(self):
        """Test that assay size patterns (without thresholds) are present."""
        result = get_all_comments()
        self.assertIn("Assay size <", result)
        self.assertIn("Assay size >", result)


class TestExplodeAssayComparability(unittest.TestCase):
    """Tests for explode_assay_comparability function."""

    def test_explode_basic(self):
        """Test basic explosion of multi-assay data."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "repeat": [0, 1],
                "activity_id": ["ACT1|ACT2|ACT3", "ACT4|ACT5"],
                "assay_chembl_id": ["ASSAY1|ASSAY2|ASSAY3", "ASSAY4|ASSAY5"],
                "pchembl_value": ["6.0|6.5|7.0", "5.5|6.0"],
                "data_processing_comment": ["||", "|"],
                "data_dropping_comment": ["||", "|"],
                "standard_type": ["IC50|IC50|IC50", "Ki|Ki"],
                "canonical_smiles": ["CCCC|CCCC|CCCC", "CCCO|CCCO"],
            }
        )

        result = explode_assay_comparability(df)

        # MOL1 has 3 assays -> 3 choose 2 = 3 pairs
        # MOL2 has 2 assays -> 2 choose 2 = 1 pair
        # Total = 4 rows
        self.assertEqual(len(result), 4)

        # Check that _x and _y columns exist
        self.assertIn("assay_chembl_id_x", result.columns)
        self.assertIn("assay_chembl_id_y", result.columns)
        self.assertIn("pchembl_value_x", result.columns)
        self.assertIn("pchembl_value_y", result.columns)
        self.assertIn("canonical_smiles_x", result.columns)
        self.assertIn("canonical_smiles_y", result.columns)

    def test_explode_filters_same_assay(self):
        """Test that comparisons of same assay are filtered out."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "target_chembl_id": ["TARGET1"],
                "repeat": [0],
                "activity_id": ["ACT1|ACT2"],
                "assay_chembl_id": ["ASSAY1|ASSAY2"],
                "pchembl_value": ["6.0|6.5"],
                "data_processing_comment": ["|"],
                "data_dropping_comment": ["|"],
                "standard_type": ["IC50|IC50"],
                "canonical_smiles": ["CCCC|CCCC"],
            }
        )

        result = explode_assay_comparability(df)

        # Should only have pairs where assay_x != assay_y
        self.assertTrue(
            (result["assay_chembl_id_x"] != result["assay_chembl_id_y"]).all(),
            "All pairs should have different assays",
        )

    def test_explode_combines_comments(self):
        """Test that processing and dropping comments are combined correctly."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "target_chembl_id": ["TARGET1"],
                "repeat": [0],
                "activity_id": ["ACT1|ACT2"],
                "assay_chembl_id": ["ASSAY1|ASSAY2"],
                "pchembl_value": ["6.0|6.5"],
                "data_processing_comment": ["Calculated pChEMBL|Salt/solvent removed"],
                "data_dropping_comment": ["Assay size < 20|"],
                "standard_type": ["IC50|IC50"],
                "canonical_smiles": ["CCCC|CCCC"],
            }
        )

        result = explode_assay_comparability(df)

        # Check that combined comments contain both original comments
        self.assertIn("processing_comment", result.columns)
        self.assertIn("dropping_comment", result.columns)

        # The first row should have both comments combined
        proc_comment = result.iloc[0]["processing_comment"]
        self.assertIn("Calculated pChEMBL", proc_comment)
        self.assertIn("Salt/solvent removed", proc_comment)

    def test_explode_empty_dataframe(self):
        """Test that empty DataFrame is handled gracefully."""
        df = pd.DataFrame(
            {
                "connectivity": [],
                "target_chembl_id": [],
                "repeat": [],
                "activity_id": [],
                "assay_chembl_id": [],
                "pchembl_value": [],
                "data_processing_comment": [],
                "data_dropping_comment": [],
                "standard_type": [],
                "canonical_smiles": [],
            }
        )

        result = explode_assay_comparability(df)

        # Should return empty DataFrame with correct columns
        self.assertEqual(len(result), 0)
        self.assertIn("assay_chembl_id_x", result.columns)
        self.assertIn("assay_chembl_id_y", result.columns)

    def test_explode_custom_separator(self):
        """Test that custom separator works correctly."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "target_chembl_id": ["TARGET1"],
                "repeat": [0],
                "activity_id": ["ACT1;ACT2"],
                "assay_chembl_id": ["ASSAY1;ASSAY2"],
                "pchembl_value": ["6.0;6.5"],
                "data_processing_comment": [";"],
                "data_dropping_comment": [";"],
                "standard_type": ["IC50;IC50"],
                "canonical_smiles": ["CCCC;CCCC"],
            }
        )

        result = explode_assay_comparability(df, sep_str=";")

        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["assay_chembl_id_x"], "ASSAY1")
        self.assertEqual(result.iloc[0]["assay_chembl_id_y"], "ASSAY2")

    def test_explode_with_single_value_columns(self):
        """Test explosion when some columns are single-valued (e.g., from --id-columns).

        When data is aggregated with --id-columns, those columns become single-valued
        rather than pipe-delimited. The function should auto-detect this and handle
        them as single-value columns.
        """
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "repeat": [0, 1],
                "activity_id": ["ACT1|ACT2|ACT3", "ACT4|ACT5"],
                "assay_chembl_id": ["ASSAY1|ASSAY2|ASSAY3", "ASSAY4|ASSAY5"],
                "pchembl_value": ["6.0|6.5|7.0", "5.5|6.0"],
                "data_processing_comment": ["||", "|"],
                "data_dropping_comment": ["||", "|"],
                # standard_type is single-valued (as if aggregated with --id-columns standard_type)
                "standard_type": ["IC50", "Ki"],
                "canonical_smiles": ["CCCC|CCCC|CCCC", "CCCO|CCCO"],
            }
        )

        result = explode_assay_comparability(df)

        # MOL1 has 3 assays -> 3 choose 2 = 3 pairs
        # MOL2 has 2 assays -> 2 choose 2 = 1 pair
        # Total = 4 rows
        self.assertEqual(len(result), 4)

        # Check that standard_type is preserved as single value (not split)
        mol1_rows = result[result["connectivity"] == "MOL1"]
        mol2_rows = result[result["connectivity"] == "MOL2"]

        # standard_type should be in columns as _x suffix since it's treated as single-value
        self.assertIn("standard_type", result.columns)
        self.assertTrue((mol1_rows["standard_type"] == "IC50").all())
        self.assertTrue((mol2_rows["standard_type"] == "Ki").all())

    def test_explode_with_custom_value_column(self):
        """Test explosion with standard_value instead of pchembl_value."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "repeat": [0, 1],
                "activity_id": ["ACT1|ACT2", "ACT3|ACT4"],
                "assay_chembl_id": ["ASSAY1|ASSAY2", "ASSAY3|ASSAY4"],
                "standard_value": ["100.0|200.0", "50.0|75.0"],
                "data_processing_comment": ["|", "|"],
                "data_dropping_comment": ["|", "|"],
                "standard_type": ["Pc|Pc", "Pc|Pc"],
                "canonical_smiles": ["CCCC|CCCC", "CCCO|CCCO"],
            }
        )

        result = explode_assay_comparability(df, value_column="standard_value")

        # MOL1 has 2 assays -> 1 pair, MOL2 has 2 assays -> 1 pair
        # Total = 2 rows
        self.assertEqual(len(result), 2)

        # Check that standard_value_x and standard_value_y columns exist
        self.assertIn("standard_value_x", result.columns)
        self.assertIn("standard_value_y", result.columns)
        # pchembl columns should NOT exist
        self.assertNotIn("pchembl_value_x", result.columns)
        self.assertNotIn("pchembl_value_y", result.columns)

        # Verify values
        mol1_row = result[result["connectivity"] == "MOL1"].iloc[0]
        self.assertEqual(mol1_row["standard_value_x"], "100.0")
        self.assertEqual(mol1_row["standard_value_y"], "200.0")


class TestDeaggregateData(unittest.TestCase):
    """Tests for deaggregate_data function."""

    def test_deaggregate_basic(self):
        """Test basic de-aggregation of pipe-delimited data."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2", "MOL3"],
                "pchembl_value": ["6.0|6.5|7.0", "5.5", "8.0|8.5"],
                "activity_id": ["ACT1|ACT2|ACT3", "ACT4", "ACT5|ACT6"],
                "assay_chembl_id": ["ASSAY1|ASSAY2|ASSAY3", "ASSAY4", "ASSAY5|ASSAY6"],
            }
        )

        result = deaggregate_data(df)

        # MOL1 should split into 3 rows, MOL2 stays 1 row, MOL3 splits into 2 rows
        # Total: 3 + 1 + 2 = 6 rows
        self.assertEqual(len(result), 6)

        # Check that MOL1 data was split correctly
        mol1_data = result[result["connectivity"] == "MOL1"].sort_values("pchembl_value")
        self.assertEqual(len(mol1_data), 3)
        self.assertEqual(mol1_data["pchembl_value"].tolist(), ["6.0", "6.5", "7.0"])
        self.assertEqual(mol1_data["activity_id"].tolist(), ["ACT1", "ACT2", "ACT3"])

        # Check that MOL2 data remained unchanged
        mol2_data = result[result["connectivity"] == "MOL2"]
        self.assertEqual(len(mol2_data), 1)
        self.assertEqual(mol2_data["pchembl_value"].values[0], "5.5")

    def test_deaggregate_no_pipes(self):
        """Test that data without pipes is returned unchanged."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "pchembl_value": ["6.0", "5.5"],
                "activity_id": ["ACT1", "ACT2"],
            }
        )

        result = deaggregate_data(df)

        # Should return same data
        self.assertEqual(len(result), 2)
        pd.testing.assert_frame_equal(result.reset_index(drop=True), df)

    def test_deaggregate_empty_dataframe(self):
        """Test that empty DataFrame is handled gracefully."""
        df = pd.DataFrame({"connectivity": [], "pchembl_value": [], "activity_id": []})

        result = deaggregate_data(df)

        # Should return empty DataFrame with same structure
        self.assertEqual(len(result), 0)
        self.assertListEqual(list(result.columns), list(df.columns))

    def test_deaggregate_mixed_columns(self):
        """Test de-aggregation with some columns having pipes and others not."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],  # No pipes
                "target_chembl_id": ["TARGET1", "TARGET2"],  # No pipes
                "pchembl_value": ["6.0|6.5", "7.0|7.5|8.0"],  # Has pipes
                "activity_id": ["ACT1|ACT2", "ACT3|ACT4|ACT5"],  # Has pipes
            }
        )

        result = deaggregate_data(df)

        # MOL1: 2 rows, MOL2: 3 rows -> Total: 5 rows
        self.assertEqual(len(result), 5)

        # Check that single-value columns are preserved
        mol1_data = result[result["connectivity"] == "MOL1"]
        self.assertTrue((mol1_data["connectivity"] == "MOL1").all())
        self.assertTrue((mol1_data["target_chembl_id"] == "TARGET1").all())

    def test_deaggregate_with_nan_values(self):
        """Test de-aggregation handles NaN values correctly."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "pchembl_value": ["6.0|6.5", None],
                "activity_id": ["ACT1|ACT2", "ACT3"],
            }
        )

        result = deaggregate_data(df)

        # MOL1 should split into 2 rows, MOL2 stays as 1 (no pipe in activity_id)
        # Total: 2 + 1 = 3 rows
        self.assertEqual(len(result), 3)

        mol1_data = result[result["connectivity"] == "MOL1"]
        self.assertEqual(len(mol1_data), 2)

    def test_deaggregate_custom_separator(self):
        """Test de-aggregation with custom separator."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "pchembl_value": ["6.0;6.5", "7.0"],
                "activity_id": ["ACT1;ACT2", "ACT3"],
            }
        )

        result = deaggregate_data(df, sep_str=";")

        # MOL1 should split into 2 rows, MOL2 stays as 1
        self.assertEqual(len(result), 3)

        mol1_data = result[result["connectivity"] == "MOL1"].sort_values("pchembl_value")
        self.assertEqual(mol1_data["pchembl_value"].tolist(), ["6.0", "6.5"])
        self.assertEqual(mol1_data["activity_id"].tolist(), ["ACT1", "ACT2"])

    def test_deaggregate_preserves_data_types(self):
        """Test that de-aggregation preserves column order and structure."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "target_chembl_id": ["TARGET1"],
                "pchembl_value": ["6.0|6.5"],
                "activity_id": ["ACT1|ACT2"],
            }
        )

        result = deaggregate_data(df)

        # Check that column order is preserved
        self.assertListEqual(list(result.columns), list(df.columns))

        # Check that we have the expected rows
        self.assertEqual(len(result), 2)

    def test_deaggregate_all_rows_aggregated(self):
        """Test de-aggregation when all rows contain pipes."""
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "pchembl_value": ["6.0|6.5", "7.0|7.5"],
                "activity_id": ["ACT1|ACT2", "ACT3|ACT4"],
            }
        )

        result = deaggregate_data(df)

        # Both rows split into 2 each -> Total: 4 rows
        self.assertEqual(len(result), 4)

        # Verify all molecules are present
        self.assertEqual(len(result[result["connectivity"] == "MOL1"]), 2)
        self.assertEqual(len(result[result["connectivity"] == "MOL2"]), 2)


class TestEnumValues(unittest.TestCase):
    """Tests for enum value correctness."""

    def test_processing_comment_enum_values(self):
        """Test that ProcessingComment enum has expected values."""
        self.assertEqual(ProcessingComment.CALCULATED_PCHEMBL.value, "Calculated pChEMBL")
        self.assertEqual(ProcessingComment.SALT_SOLVENT_REMOVED.value, "Salt/solvent removed")
        self.assertEqual(
            ProcessingComment.PCHEMBL_DUPLICATION_ACROSS_DOCUMENTS.value,
            "pChEMBL Duplication Across Documents",
        )

    def test_dropping_comment_enum_values(self):
        """Test that DroppingComment enum has expected values."""
        self.assertEqual(DroppingComment.DATA_VALIDITY_COMMENT.value, "Data Validity Comment Present")
        self.assertEqual(DroppingComment.POTENTIAL_DUPLICATE.value, "Potential Duplicate")
        self.assertEqual(DroppingComment.UNDEFINED_STEREOCHEMISTRY.value, "Undefined Stereochemistry")
        self.assertEqual(DroppingComment.MUTATION_KEYWORD.value, "Mutation keyword in assay description")
        self.assertEqual(DroppingComment.UNIT_ANNOTATION_ERROR.value, "Unit Annotation Error")

    def test_assay_size_enum_values_are_patterns(self):
        """Test that assay size enum values are patterns (without thresholds)."""
        self.assertEqual(DroppingComment.ASSAY_SIZE_TOO_SMALL.value, "Assay size <")
        self.assertEqual(DroppingComment.ASSAY_SIZE_TOO_LARGE.value, "Assay size >")

    def test_enum_string_inheritance(self):
        """Test that enums inherit from str for easy comparison."""
        self.assertIsInstance(ProcessingComment.CALCULATED_PCHEMBL, str)
        self.assertIsInstance(DroppingComment.POTENTIAL_DUPLICATE, str)


if __name__ == "__main__":
    unittest.main()
