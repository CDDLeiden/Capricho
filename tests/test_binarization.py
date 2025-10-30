import unittest

import numpy as np
import pandas as pd

from Capricho.core.binarization import (
    invert_relation_for_pchembl,
    binarize_aggregated_data,
)
from Capricho.core.default_fields import DATA_DROPPING_COMMENT


class TestInvertRelationForPchembl(unittest.TestCase):
    def test_invert_less_than(self):
        self.assertEqual(invert_relation_for_pchembl("<"), ">")

    def test_invert_greater_than(self):
        self.assertEqual(invert_relation_for_pchembl(">"), "<")

    def test_invert_equals(self):
        self.assertEqual(invert_relation_for_pchembl("="), "=")

    def test_invert_less_than_equals(self):
        self.assertEqual(invert_relation_for_pchembl("<="), ">=")

    def test_invert_greater_than_equals(self):
        self.assertEqual(invert_relation_for_pchembl(">="), "<=")

    def test_invalid_relation(self):
        with self.assertRaises(ValueError):
            invert_relation_for_pchembl("invalid")

    def test_invert_approximately_equal(self):
        self.assertEqual(invert_relation_for_pchembl("~"), "~")

    def test_invert_much_greater_than(self):
        self.assertEqual(invert_relation_for_pchembl(">>"), "<<")

    def test_invert_much_less_than(self):
        self.assertEqual(invert_relation_for_pchembl("<<"), ">>")


class TestBinarizeAggregatedData(unittest.TestCase):
    def setUp(self):
        self.threshold = 6.0

    def test_basic_discrete_binarization(self):
        """Test basic binarization with only '=' relations"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 5.0, 6.0],
                "standard_relation": ["=", "=", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # Check binary labels
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # 7.0 >= 6.0
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # 5.0 < 6.0
        self.assertEqual(result.loc[2, "activity_binary"], 1)  # 6.0 >= 6.0 (edge case)

    def test_censored_active_binarization(self):
        """Test binarization with '<' (censored active) relations"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 5.0, 6.0],
                "standard_relation": ["<", "<", "<"],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # '<' means compound is MORE active (lower concentration)
        # If pchembl >= threshold, definitely active
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # 7.0 >= 6.0 → active
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # 5.0 < 6.0 → might be inactive
        self.assertEqual(result.loc[2, "activity_binary"], 1)  # 6.0 >= 6.0 → active

    def test_censored_inactive_binarization(self):
        """Test binarization with '>' (censored inactive) relations"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 5.0, 6.0],
                "standard_relation": [">", ">", ">"],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # '>' means compound is LESS active (higher concentration)
        # If pchembl <= threshold, definitely inactive
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # 7.0 > 6.0 → might be active
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # 5.0 <= 6.0 → inactive
        self.assertEqual(result.loc[2, "activity_binary"], 0)  # 6.0 <= 6.0 → inactive

    def test_mixed_relations_with_agreement(self):
        """Test mixed discrete and censored measurements that agree"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 8.0],
                "standard_relation": ["<", "="],  # censored active + discrete
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # Both should be active
        self.assertEqual(result.loc[0, "activity_binary"], 1)
        self.assertEqual(result.loc[1, "activity_binary"], 1)

        # Should not have conflict flag
        if DATA_DROPPING_COMMENT in result.columns:
            self.assertTrue(
                result[DATA_DROPPING_COMMENT].isna().all()
                or not result[DATA_DROPPING_COMMENT].str.contains("Non-agreeing").any()
            )

    def test_mixed_relations_with_disagreement(self):
        """Test mixed discrete and censored measurements that disagree"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 5.0],
                "standard_relation": ["<", "="],  # censored says active, discrete says inactive
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # Should have conflict flag for both rows
        self.assertIn(DATA_DROPPING_COMMENT, result.columns)
        self.assertTrue(result[DATA_DROPPING_COMMENT].str.contains("Non-agreeing", na=False).all())

    def test_nan_values(self):
        """Test handling of NaN values in pchembl column"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, np.nan, 5.0],
                "standard_relation": ["=", "=", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # Check that NaN is preserved
        self.assertEqual(result.loc[0, "activity_binary"], 1)
        self.assertTrue(pd.isna(result.loc[1, "activity_binary"]))
        self.assertEqual(result.loc[2, "activity_binary"], 0)

    def test_missing_standard_relation_column(self):
        """Test that function handles missing standard_relation column"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 5.0],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # Should assume '=' for all rows
        self.assertEqual(result.loc[0, "activity_binary"], 1)
        self.assertEqual(result.loc[1, "activity_binary"], 0)

    def test_custom_column_names(self):
        """Test using custom column names"""
        df = pd.DataFrame(
            {
                "my_compound_id": ["CONN1", "CONN2"],
                "my_target_id": ["TARGET1", "TARGET1"],
                "my_value": [7.0, 5.0],
                "my_relation": ["=", "="],
            }
        )

        result = binarize_aggregated_data(
            df,
            threshold=self.threshold,
            value_column="my_value",
            compound_id_col="my_compound_id",
            target_id_col="my_target_id",
            relation_col="my_relation",
            output_binary_col="my_binary",
        )

        self.assertIn("my_binary", result.columns)
        self.assertEqual(result.loc[0, "my_binary"], 1)
        self.assertEqual(result.loc[1, "my_binary"], 0)

    def test_different_threshold(self):
        """Test with non-default threshold"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 6.5, 6.0],
                "standard_relation": ["=", "=", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=6.5)

        self.assertEqual(result.loc[0, "activity_binary"], 1)  # 7.0 >= 6.5
        self.assertEqual(result.loc[1, "activity_binary"], 1)  # 6.5 >= 6.5
        self.assertEqual(result.loc[2, "activity_binary"], 0)  # 6.0 < 6.5

    def test_median_column(self):
        """Test using median instead of mean"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 5.0],
                "pchembl_value_median": [6.5, 5.5],
                "standard_relation": ["=", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold, value_column="pchembl_value_median")

        # Using median values now
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # 6.5 >= 6.0
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # 5.5 < 6.0

    def test_multiple_compound_target_pairs(self):
        """Test with multiple compound-target pairs"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1", "CONN2", "CONN2"],
                "target_chembl_id": ["TARGET1", "TARGET2", "TARGET1", "TARGET2"],
                "pchembl_value_mean": [7.0, 5.0, 6.5, 5.5],
                "standard_relation": ["=", "=", "=", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # Each compound-target pair should be binarized independently
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # CONN1-TARGET1
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # CONN1-TARGET2
        self.assertEqual(result.loc[2, "activity_binary"], 1)  # CONN2-TARGET1
        self.assertEqual(result.loc[3, "activity_binary"], 0)  # CONN2-TARGET2

    def test_complex_mixed_relations(self):
        """Test complex scenario with multiple relations for same compound-target"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1", "CONN1"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [5.5, 6.5, 7.5],
                "standard_relation": ["<", "=", ">"],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # < at 5.5: might be inactive (below threshold)
        self.assertEqual(result.loc[0, "activity_binary"], 0)
        # = at 6.5: active (above threshold)
        self.assertEqual(result.loc[1, "activity_binary"], 1)
        # > at 7.5: might be active (above threshold)
        self.assertEqual(result.loc[2, "activity_binary"], 1)

    def test_approximately_equal_binarization(self):
        """Test binarization with '~' (approximately equal) relation"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3", "CONN4"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 6.5, 6.0, 5.5],
                "standard_relation": ["~", "~", "~", "~"],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # ~ at 7.0: lower bound = 6.5 >= 6.0 → active
        self.assertEqual(result.loc[0, "activity_binary"], 1)
        # ~ at 6.5: lower bound = 6.0 >= 6.0 → active (edge case)
        self.assertEqual(result.loc[1, "activity_binary"], 1)
        # ~ at 6.0: lower bound = 5.5 < 6.0 → inactive
        self.assertEqual(result.loc[2, "activity_binary"], 0)
        # ~ at 5.5: lower bound = 5.0 < 6.0 → inactive
        self.assertEqual(result.loc[3, "activity_binary"], 0)

    def test_much_greater_than_binarization(self):
        """Test binarization with '>>' (much greater than) relation"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 6.0, 5.0],
                "standard_relation": [">>", ">>", ">>"],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # >> means compound is MUCH LESS active (concentration >> reported)
        # Treat like '>': if pchembl <= threshold → inactive
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # 7.0 > 6.0 → might be active
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # 6.0 <= 6.0 → inactive
        self.assertEqual(result.loc[2, "activity_binary"], 0)  # 5.0 < 6.0 → inactive

    def test_much_less_than_binarization(self):
        """Test binarization with '<<' (much less than) relation"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN2", "CONN3"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "pchembl_value_mean": [7.0, 6.0, 5.0],
                "standard_relation": ["<<", "<<", "<<"],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # << means compound is MUCH MORE active (concentration << reported)
        # Treat like '<': if pchembl >= threshold → active
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # 7.0 >= 6.0 → active
        self.assertEqual(result.loc[1, "activity_binary"], 1)  # 6.0 >= 6.0 → active (edge case)
        self.assertEqual(result.loc[2, "activity_binary"], 0)  # 5.0 < 6.0 → might be inactive

    def test_mixed_with_approximate(self):
        """Test mixed '~' and '=' measurements for same compound-target"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "pchembl_value_mean": [6.8, 7.0],
                "standard_relation": ["~", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold)

        # Both should be active
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # ~ at 6.8: lower bound 6.3
        self.assertEqual(result.loc[1, "activity_binary"], 1)  # = at 7.0

    def test_different_mutants_not_flagged_when_compare_across_mutants_false(self):
        """Test that different mutants are NOT considered as disagreements when compare_across_mutants=False"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "mutation": ["Wild-type", "V600E"],
                "pchembl_value_mean": [7.0, 5.0],
                "standard_relation": ["=", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold, compare_across_mutants=False)

        # Both should be binarized
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # Active on wild-type
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # Inactive on mutant

        # Should NOT have conflict flag because they are different mutants
        if DATA_DROPPING_COMMENT in result.columns:
            self.assertTrue(
                result[DATA_DROPPING_COMMENT].isna().all()
                or not result[DATA_DROPPING_COMMENT].str.contains("Non-agreeing", na=False).any()
            )

    def test_different_mutants_flagged_when_compare_across_mutants_true(self):
        """Test that different mutants ARE considered as disagreements when compare_across_mutants=True"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "mutation": ["Wild-type", "V600E"],
                "pchembl_value_mean": [7.0, 5.0],
                "standard_relation": ["=", "="],
            }
        )

        result = binarize_aggregated_data(df, threshold=self.threshold, compare_across_mutants=True)

        # Both should be binarized
        self.assertEqual(result.loc[0, "activity_binary"], 1)  # Active on wild-type
        self.assertEqual(result.loc[1, "activity_binary"], 0)  # Inactive on mutant

        # SHOULD have conflict flag because we're aggregating across mutants
        self.assertIn(DATA_DROPPING_COMMENT, result.columns)
        self.assertTrue(result[DATA_DROPPING_COMMENT].str.contains("Non-agreeing", na=False).all())

    def test_same_mutant_with_disagreement_always_flagged(self):
        """Test that disagreements on the same mutant are always flagged regardless of compare_across_mutants"""
        df = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1"],
                "target_chembl_id": ["TARGET1", "TARGET1"],
                "mutation": ["Wild-type", "Wild-type"],
                "pchembl_value_mean": [7.0, 5.0],
                "standard_relation": ["<", "="],  # Censored says active, discrete says inactive
            }
        )

        # Test with compare_across_mutants=False
        result_false = binarize_aggregated_data(df, threshold=self.threshold, compare_across_mutants=False)
        self.assertIn(DATA_DROPPING_COMMENT, result_false.columns)
        self.assertTrue(result_false[DATA_DROPPING_COMMENT].str.contains("Non-agreeing", na=False).all())

        # Test with compare_across_mutants=True (should have same result)
        result_true = binarize_aggregated_data(df, threshold=self.threshold, compare_across_mutants=True)
        self.assertIn(DATA_DROPPING_COMMENT, result_true.columns)
        self.assertTrue(result_true[DATA_DROPPING_COMMENT].str.contains("Non-agreeing", na=False).all())


if __name__ == "__main__":
    unittest.main()
