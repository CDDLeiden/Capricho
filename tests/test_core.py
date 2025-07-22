import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from Capricho.core import pandas_helper, smiles_utils, stats_make
from Capricho.core.pandas_helper import conflicting_duplicates


class TestPandasHelper(unittest.TestCase):
    def test_format_value(self):
        self.assertEqual(pandas_helper.format_value(3.14159), "3.14")
        self.assertEqual(pandas_helper.format_value(2000), "2000")
        self.assertEqual(pandas_helper.format_value("test"), "test")

    def test_aggr_val_series(self):
        series = pd.Series([1, 2.5, 3, 4.7])
        self.assertEqual(pandas_helper.aggr_val_series(series), "1.00|2.50|3.00|4.70")

    def test_get_mad(self):
        values = [1, 2, 3, 4, 5]
        self.assertEqual(pandas_helper.get_mad(values), 1.0)
        self.assertTrue(np.isnan(pandas_helper.get_mad([1])))

    def test_apply_func_grpd(self):
        df = pd.DataFrame({"A": [1, 1, 2, 2], "B": [1, 2, 3, 4]})
        grouped = df.groupby("A")
        result = pandas_helper.apply_func_grpd(grouped, pandas_helper.aggr_val_series, ["A"], "B")
        expected = pd.DataFrame({"A": [1, 2], "B": ["1|2", "3|4"]})
        pd.testing.assert_frame_equal(result, expected)

    def test_assign_stats(self):
        df = pd.DataFrame({"value": ["1|2|3", "4|5|6"]})
        result = pandas_helper.assign_stats(df, sep="|", value_col="value")
        self.assertIn("value_mean", result.columns)
        self.assertIn("value_std", result.columns)
        self.assertIn("value_median", result.columns)
        self.assertIn("value_counts", result.columns)

        self.assertEqual(result["value_mean"].iloc[0], 2.0)
        self.assertAlmostEquals(result["value_std"].iloc[0], 0.816497, 6)
        self.assertEqual(result["value_median"].iloc[0], 2.0)
        self.assertEqual(result["value_counts"].iloc[0], 3)

        self.assertEqual(result["value_mean"].iloc[1], 5.0)
        self.assertAlmostEquals(result["value_std"].iloc[1], 0.816497, 6)
        self.assertEqual(result["value_median"].iloc[1], 5.0)
        self.assertEqual(result["value_counts"].iloc[1], 3)

    def test_conflicting_duplicates(self):
        data = {
            "A": [
                1,
                1,
                1,
                2,
                2,
            ],
            "B": [
                "x",
                "x",
                "x",
                "y",
                "y",
            ],
            "C": [
                "p",
                "p",
                "p",
                "q",
                "q",
            ],
            "D": [
                10,  # Same A, B, C (e.g. mol identifiers) but different D (document)
                10,  # Same as above
                20,  # Same as above
                30,  # Same everything -> reported in same document; shouldn't drop
                30,  # Same as above
            ],
            "year": [
                2020,
                2020,
                2021,
                2021,
                2021,
            ],
        }
        df = pd.DataFrame(data).sort_values(by="year")
        mask = conflicting_duplicates(df, key_subset=["A", "B", "C"], diff_subset=["D"])
        expected_flags = [True, True, True, False, False]
        true_idxs = np.where(mask)[0]
        assert np.array_equal(
            true_idxs, np.array([0, 1, 2])
        ), "Expected indices with conflicting duplicates do not match."
        pd.testing.assert_series_equal(mask, pd.Series(expected_flags, index=df.index))


class TestSmilesUtils(unittest.TestCase):
    def test_clean_mixtures(self):
        self.assertEqual(smiles_utils.clean_mixtures("CC.Cl"), "CC")
        self.assertEqual(smiles_utils.clean_mixtures("CC.Na+"), "CC")
        self.assertEqual(smiles_utils.clean_mixtures("CC.O"), "CC")


class TestStatsMake(unittest.TestCase):
    def setUp(self):
        self.testroot = Path(__file__).parent
        self.not_aggr_df = pd.read_csv(self.testroot / "resources/ADORA3_data_not_aggregated.csv")
        self.aggr_df = pd.read_csv(self.testroot / "resources/ADORA3_data.csv")

    # note - test_process_repeat_mols not tested here but indirectly on `test_workflow.py`

    def test_repeated_indices_from_array_series(self):
        series = pd.Series([np.array([1, 2]), np.array([1, 2]), np.array([3, 4])])
        result = stats_make.repeated_indices_from_array_series(series)
        self.assertEqual(result, [[0, 1]])


if __name__ == "__main__":
    unittest.main()
