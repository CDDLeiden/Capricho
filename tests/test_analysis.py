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


class TestDeduplicateAggregatedValues(unittest.TestCase):
    """Tests for deduplicate_aggregated_values function."""

    def test_deduplicate_identical_values(self):
        """Test that identical pchembl values within a row are deduplicated."""
        from Capricho.analysis import deduplicate_aggregated_values

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "pchembl_value": ["8.00|8.00|8.00|6.92"],
                "activity_id": ["ACT1|ACT2|ACT3|ACT4"],
                "document_chembl_id": ["DOC1|DOC2|DOC3|DOC4"],
                "data_processing_comment": [
                    "pChEMBL Duplication Across Documents|pChEMBL Duplication Across Documents||"
                ],
            }
        )

        result = deduplicate_aggregated_values(df)

        # Should keep only unique pchembl values: 8.00 and 6.92
        self.assertEqual(len(result), 1)
        pchembl_vals = result.iloc[0]["pchembl_value"].split("|")
        self.assertEqual(len(pchembl_vals), 2)
        self.assertIn("8.00", pchembl_vals)
        self.assertIn("6.92", pchembl_vals)

    def test_deduplicate_keeps_first_occurrence(self):
        """Test that deduplication keeps first occurrence of each value."""
        from Capricho.analysis import deduplicate_aggregated_values

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "pchembl_value": ["8.00|6.50|8.00"],
                "activity_id": ["ACT1|ACT2|ACT3"],
                "document_chembl_id": ["DOC1|DOC2|DOC3"],
                "data_processing_comment": ["||"],
            }
        )

        result = deduplicate_aggregated_values(df)

        # First 8.00 should be kept (ACT1, DOC1)
        pchembl_vals = result.iloc[0]["pchembl_value"].split("|")
        activity_ids = result.iloc[0]["activity_id"].split("|")
        doc_ids = result.iloc[0]["document_chembl_id"].split("|")

        self.assertEqual(len(pchembl_vals), 2)
        # Check order: 8.00 first, then 6.50
        self.assertEqual(pchembl_vals[0], "8.00")
        self.assertEqual(pchembl_vals[1], "6.50")
        self.assertEqual(activity_ids[0], "ACT1")
        self.assertEqual(doc_ids[0], "DOC1")

    def test_deduplicate_no_duplicates(self):
        """Test that rows without duplicates are unchanged."""
        from Capricho.analysis import deduplicate_aggregated_values

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "pchembl_value": ["8.00|7.50|6.92"],
                "activity_id": ["ACT1|ACT2|ACT3"],
                "document_chembl_id": ["DOC1|DOC2|DOC3"],
                "data_processing_comment": ["||"],
            }
        )

        result = deduplicate_aggregated_values(df)

        # No duplicates, should remain unchanged
        self.assertEqual(result.iloc[0]["pchembl_value"], "8.00|7.50|6.92")
        self.assertEqual(result.iloc[0]["activity_id"], "ACT1|ACT2|ACT3")

    def test_deduplicate_single_value_rows(self):
        """Test that single-value rows are unchanged."""
        from Capricho.analysis import deduplicate_aggregated_values

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "pchembl_value": ["8.00", "7.50|7.50"],
                "activity_id": ["ACT1", "ACT2|ACT3"],
                "document_chembl_id": ["DOC1", "DOC2|DOC3"],
                "data_processing_comment": ["", "|"],
            }
        )

        result = deduplicate_aggregated_values(df)

        # MOL1 unchanged (single value)
        mol1 = result[result["connectivity"] == "MOL1"].iloc[0]
        self.assertEqual(mol1["pchembl_value"], "8.00")

        # MOL2 deduplicated (7.50|7.50 -> 7.50)
        mol2 = result[result["connectivity"] == "MOL2"].iloc[0]
        self.assertEqual(mol2["pchembl_value"], "7.50")

    def test_deduplicate_preserves_other_columns(self):
        """Test that columns not involved in deduplication are preserved."""
        from Capricho.analysis import deduplicate_aggregated_values

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "target_chembl_id": ["TARGET1"],
                "smiles": ["CCCC"],
                "pchembl_value": ["8.00|8.00"],
                "activity_id": ["ACT1|ACT2"],
                "document_chembl_id": ["DOC1|DOC2"],
                "data_processing_comment": ["|"],
                "pchembl_value_mean": [8.0],
                "pchembl_value_counts": [2.0],
            }
        )

        result = deduplicate_aggregated_values(df)

        # Single-value columns should be preserved
        self.assertEqual(result.iloc[0]["connectivity"], "MOL1")
        self.assertEqual(result.iloc[0]["target_chembl_id"], "TARGET1")
        self.assertEqual(result.iloc[0]["smiles"], "CCCC")

    def test_deduplicate_empty_dataframe(self):
        """Test that empty DataFrame is handled gracefully."""
        from Capricho.analysis import deduplicate_aggregated_values

        df = pd.DataFrame(
            {
                "connectivity": [],
                "pchembl_value": [],
                "activity_id": [],
                "document_chembl_id": [],
                "data_processing_comment": [],
            }
        )

        result = deduplicate_aggregated_values(df)

        self.assertEqual(len(result), 0)


class TestRecalculateAggregatedStats(unittest.TestCase):
    """Tests for recalculate_aggregated_stats function."""

    def test_recalculate_stats_after_dedup(self):
        """Test that stats are recalculated after deduplication."""
        from Capricho.analysis import recalculate_aggregated_stats

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "pchembl_value": ["8.00|6.00"],  # Already deduplicated
                "pchembl_value_mean": [7.33],  # Old mean from 8.00|8.00|6.00
                "pchembl_value_median": [8.0],
                "pchembl_value_std": [1.15],
                "pchembl_value_counts": [3.0],
            }
        )

        result = recalculate_aggregated_stats(df)

        # New mean should be (8.0 + 6.0) / 2 = 7.0
        self.assertAlmostEqual(result.iloc[0]["pchembl_value_mean"], 7.0, places=2)
        # New median should be 7.0
        self.assertAlmostEqual(result.iloc[0]["pchembl_value_median"], 7.0, places=2)
        # New count should be 2
        self.assertEqual(result.iloc[0]["pchembl_value_counts"], 2)

    def test_recalculate_stats_single_value(self):
        """Test that single values get correct stats."""
        from Capricho.analysis import recalculate_aggregated_stats

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1"],
                "pchembl_value": ["8.00"],
                "pchembl_value_mean": [8.0],
                "pchembl_value_median": [8.0],
                "pchembl_value_std": [0.0],
                "pchembl_value_counts": [1.0],
            }
        )

        result = recalculate_aggregated_stats(df)

        self.assertAlmostEqual(result.iloc[0]["pchembl_value_mean"], 8.0)
        self.assertAlmostEqual(result.iloc[0]["pchembl_value_median"], 8.0)
        self.assertEqual(result.iloc[0]["pchembl_value_counts"], 1)

    def test_recalculate_stats_multiple_rows(self):
        """Test recalculation works for multiple rows."""
        from Capricho.analysis import recalculate_aggregated_stats

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL2"],
                "pchembl_value": ["8.00|6.00|7.00", "5.50"],
                "pchembl_value_mean": [7.0, 5.5],
                "pchembl_value_median": [7.0, 5.5],
                "pchembl_value_std": [1.0, 0.0],
                "pchembl_value_counts": [3.0, 1.0],
            }
        )

        result = recalculate_aggregated_stats(df)

        # MOL1: mean of 8, 6, 7 = 7.0
        self.assertAlmostEqual(result.iloc[0]["pchembl_value_mean"], 7.0, places=2)
        self.assertEqual(result.iloc[0]["pchembl_value_counts"], 3)

        # MOL2: single value, unchanged
        self.assertAlmostEqual(result.iloc[1]["pchembl_value_mean"], 5.5)
        self.assertEqual(result.iloc[1]["pchembl_value_counts"], 1)

    def test_recalculate_stats_empty_dataframe(self):
        """Test that empty DataFrame is handled gracefully."""
        from Capricho.analysis import recalculate_aggregated_stats

        df = pd.DataFrame(
            {
                "connectivity": [],
                "pchembl_value": [],
                "pchembl_value_mean": [],
                "pchembl_value_median": [],
                "pchembl_value_std": [],
                "pchembl_value_counts": [],
            }
        )

        result = recalculate_aggregated_stats(df)

        self.assertEqual(len(result), 0)


class TestResolveAnnotationErrors(unittest.TestCase):
    """Tests for resolve_annotation_errors function."""

    def test_resolve_keeps_earliest_document(self):
        """Test that resolution keeps measurement from earliest document."""
        from Capricho.analysis import resolve_annotation_errors

        # Two measurements for same molecule, different assays, differ by 3.0
        # Value 5.0 from 2010, value 8.0 from 2015 - should keep 5.0
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL1"],
                "pchembl_value": ["5.00", "8.00"],
                "year": ["2010", "2015"],
                "assay_chembl_id": ["ASSAY1", "ASSAY2"],
                "molecule_chembl_id": ["CHEMBL123", "CHEMBL123"],
                "activity_id": ["ACT1", "ACT2"],
            }
        )

        result = resolve_annotation_errors(df, strategy="first")

        # Should keep only the 2010 measurement (value=5.0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["pchembl_value"], "5.00")
        self.assertEqual(result.iloc[0]["year"], "2010")

    def test_resolve_6_log_unit_difference(self):
        """Test that 6.0 log unit differences are also detected."""
        from Capricho.analysis import resolve_annotation_errors

        # Differ by 6.0: 4.0 vs 10.0
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL1"],
                "pchembl_value": ["4.00", "10.00"],
                "year": ["2012", "2008"],
                "assay_chembl_id": ["ASSAY1", "ASSAY2"],
                "molecule_chembl_id": ["CHEMBL123", "CHEMBL123"],
                "activity_id": ["ACT1", "ACT2"],
            }
        )

        result = resolve_annotation_errors(df, strategy="first")

        # Should keep 2008 measurement (value=10.0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["pchembl_value"], "10.00")
        self.assertEqual(result.iloc[0]["year"], "2008")

    def test_resolve_no_error_unchanged(self):
        """Test that measurements without annotation errors are unchanged."""
        from Capricho.analysis import resolve_annotation_errors

        # Two measurements differ by 1.5 (not 3.0 or 6.0) - no error
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL1"],
                "pchembl_value": ["5.00", "6.50"],
                "year": ["2010", "2015"],
                "assay_chembl_id": ["ASSAY1", "ASSAY2"],
                "molecule_chembl_id": ["CHEMBL123", "CHEMBL123"],
                "activity_id": ["ACT1", "ACT2"],
            }
        )

        result = resolve_annotation_errors(df, strategy="first")

        # Both should be kept
        self.assertEqual(len(result), 2)

    def test_resolve_same_assay_ignored(self):
        """Test that pairs from the same assay are not flagged."""
        from Capricho.analysis import resolve_annotation_errors

        # Same assay - even with 3.0 difference, not an annotation error
        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL1"],
                "pchembl_value": ["5.00", "8.00"],
                "year": ["2010", "2015"],
                "assay_chembl_id": ["ASSAY1", "ASSAY1"],  # Same assay
                "molecule_chembl_id": ["CHEMBL123", "CHEMBL123"],
                "activity_id": ["ACT1", "ACT2"],
            }
        )

        result = resolve_annotation_errors(df, strategy="first")

        # Both should be kept (same assay doesn't trigger annotation error)
        self.assertEqual(len(result), 2)

    def test_resolve_multiple_molecules(self):
        """Test resolution works across multiple molecules."""
        from Capricho.analysis import resolve_annotation_errors

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL1", "MOL2", "MOL2"],
                "pchembl_value": ["5.00", "8.00", "6.00", "7.00"],
                "year": ["2010", "2015", "2012", "2014"],
                "assay_chembl_id": ["ASSAY1", "ASSAY2", "ASSAY3", "ASSAY4"],
                "molecule_chembl_id": ["CHEMBL1", "CHEMBL1", "CHEMBL2", "CHEMBL2"],
                "activity_id": ["ACT1", "ACT2", "ACT3", "ACT4"],
            }
        )

        result = resolve_annotation_errors(df, strategy="first")

        # MOL1 has 3.0 diff -> keep earliest (2010, value=5.0)
        # MOL2 has 1.0 diff -> keep both
        self.assertEqual(len(result), 3)
        mol1_vals = result[result["connectivity"] == "MOL1"]["pchembl_value"].tolist()
        mol2_vals = result[result["connectivity"] == "MOL2"]["pchembl_value"].tolist()
        self.assertEqual(mol1_vals, ["5.00"])
        self.assertEqual(sorted(mol2_vals), ["6.00", "7.00"])

    def test_resolve_nan_years_excluded(self):
        """Test that measurements with nan years are handled gracefully."""
        from Capricho.analysis import resolve_annotation_errors

        df = pd.DataFrame(
            {
                "connectivity": ["MOL1", "MOL1"],
                "pchembl_value": ["5.00", "8.00"],
                "year": ["nan", "2015"],
                "assay_chembl_id": ["ASSAY1", "ASSAY2"],
                "molecule_chembl_id": ["CHEMBL123", "CHEMBL123"],
                "activity_id": ["ACT1", "ACT2"],
            }
        )

        result = resolve_annotation_errors(df, strategy="first")

        # Should keep 2015 (the one with valid year) since nan is not comparable
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["year"], "2015")


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


class TestLogScaleTransformation(unittest.TestCase):
    """Tests for log scale transformation with scale factor support."""

    def test_log_transform_with_scale_factor_produces_positive_values(self):
        """Test that log transformation with scale factor produces correct positive values.

        When data is in units of 10^-6 cm/s (stored as coefficients like 5 meaning
        5 × 10^-6 cm/s), the transformation -log10(value * 1e-6) should produce
        positive values around 5-6, not negative values.
        """
        import numpy as np

        from Capricho.analysis import plot_subset

        # Typical Caco-2 permeability values in 10^-6 cm/s units
        df = pd.DataFrame(
            {
                "standard_value_x": [5.0, 10.0, 50.0],  # 5, 10, 50 × 10^-6 cm/s
                "standard_value_y": [4.0, 12.0, 45.0],
            }
        )

        # Expected transformation: -log10(5 * 1e-6) = -log10(5e-6) = 5.3
        # For value=5 with scale_factor=1e-6:
        # -log10(5) + (-log10(1e-6)) = -0.7 + 6 = 5.3
        fig, ax = plot_subset(
            df,
            value_column="standard_value",
            log_transform=True,
            log_scale_factor=1e-6,
        )

        # Get the plotted data from the scatter collection
        scatter = ax.collections[0]
        offsets = scatter.get_offsets()

        # All values should be positive (in the 4-7 range for typical Caco-2 data)
        self.assertTrue((offsets[:, 0] > 0).all(), "X values should be positive after transformation")
        self.assertTrue((offsets[:, 1] > 0).all(), "Y values should be positive after transformation")

        # Verify specific values: -log10(5e-6) ≈ 5.30
        expected_x_first = -np.log10(5.0 * 1e-6)  # ≈ 5.30
        self.assertAlmostEqual(offsets[0, 0], expected_x_first, places=2)

        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_log_transform_without_scale_factor_unchanged(self):
        """Test that log transformation without scale factor (default) works as before."""
        import numpy as np

        from Capricho.analysis import plot_subset

        df = pd.DataFrame(
            {
                "standard_value_x": [100.0, 1000.0],
                "standard_value_y": [200.0, 800.0],
            }
        )

        fig, ax = plot_subset(
            df,
            value_column="standard_value",
            log_transform=True,
            # No log_scale_factor specified - should default to 1.0
        )

        scatter = ax.collections[0]
        offsets = scatter.get_offsets()

        # With scale_factor=1.0 and epsilon=1e-9, result is essentially -log10(value)
        # For value=100: -log10(100 + 1e-9) ≈ -2.0
        expected_x_first = -np.log10(100.0 + 1e-9)
        self.assertAlmostEqual(offsets[0, 0], expected_x_first, places=2)

        import matplotlib.pyplot as plt

        plt.close(fig)


if __name__ == "__main__":
    unittest.main()
