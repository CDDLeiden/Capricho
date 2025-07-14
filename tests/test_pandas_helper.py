import unittest

import numpy as np
import pandas as pd

from CompoundMapper.core.pandas_helper import conflicting_duplicates


class TestActivityCuration(unittest.TestCase):

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


if __name__ == "__main__":
    unittest.main()
