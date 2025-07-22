"""Module to test the workflow functions within /cli. To reproduce the same data in the
`resources/` directory, check the check `ADORA3_data_recipe.json`.
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from Capricho.cli.chembl_data_pipeline import aggregate_data


class TestFetchFromChEMBL(unittest.TestCase):
    def setUp(self):
        self.testroot = Path(__file__).parent
        self.aggr_df = pd.read_csv(self.testroot / "resources/ADORA3_data.csv")
        self.not_aggr_df = pd.read_csv(self.testroot / "resources/ADORA3_data_not_aggregated.csv")

    def test_aggregate_data(self):
        aggr_df = aggregate_data(
            self.not_aggr_df,
            chirality=True,
            metadata_cols=[],
            output_path=None,
            compound_equality="connectivity",
        )
        with tempfile.TemporaryDirectory() as tmpdirname:
            outfile = Path(tmpdirname) / "ADORA3_data.csv"
            aggr_df.to_csv(outfile, index=False)
            aggr_df = pd.read_csv(outfile)

        aggr_df.sort_values(["activity_id"], inplace=True, ignore_index=True)
        self.aggr_df.sort_values(["activity_id"], inplace=True, ignore_index=True)

        # assert if the columns are the same
        np.testing.assert_array_equal(self.aggr_df.columns, aggr_df.columns)

        # Compare if the bioactivity data is the same
        np.testing.assert_array_equal(aggr_df.pchembl_value.values, self.aggr_df.pchembl_value.values)
        np.testing.assert_array_equal(
            aggr_df.pchembl_value_median.values, self.aggr_df.pchembl_value_median.values
        )
        np.testing.assert_array_equal(
            aggr_df.pchembl_value_mean.values, self.aggr_df.pchembl_value_mean.values
        )
        np.testing.assert_array_equal(
            aggr_df.molecule_chembl_id.values, self.aggr_df.molecule_chembl_id.values
        )

        # Compare if the assays are the same
        np.testing.assert_array_equal(aggr_df.assay_chembl_id.values, self.aggr_df.assay_chembl_id.values)
        np.testing.assert_array_equal(aggr_df.activity_id.values, self.aggr_df.activity_id.values)


if __name__ == "__main__":
    unittest.main()
