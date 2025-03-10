"""Module to test the workflow functions within /cli. To reproduce the same data found in test/resources, run:

`getchembl -tids CHEMBL256 -o ADORA3_data.csv -c 7 8 9 -biotype Ki -chiral -v 28 -noaggr`
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from CompoundMapper.cli.chembl_data_pipeline import aggregate_data


class TestFetchFromChEMBL(unittest.TestCase):
    def setUp(self):
        self.testroot = Path(__file__).parent
        self.aggr_df = pd.read_csv(self.testroot / "resources/ADORA3_data.csv")
        self.not_aggr_df = pd.read_csv(self.testroot / "resources/ADORA3_data_not_aggregated.csv")

    def test_aggregate_data(self):
        aggr_df = aggregate_data(
            self.not_aggr_df,
            chirality=False,
            chembl_version=31,
            metadata_cols=[],
            output_path=None,
        )
        with tempfile.TemporaryDirectory() as tmpdirname:
            outfile = Path(tmpdirname) / "ADORA3_data.csv"
            aggr_df.to_csv(outfile, index=False)
            aggr_df = pd.read_csv(outfile)

        aggr_df.sort_values(["activity_id"], inplace=True, ignore_index=True)
        self.aggr_df.sort_values(["activity_id"], inplace=True, ignore_index=True)

        # assert if the columns are the same
        np.testing.assert_array_equal(aggr_df.columns, self.aggr_df.columns)

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
