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

    def test_censored_measurement_aggregation(self):
        """Test that censored measurements are only aggregated if they have identical relation AND value."""
        test_data = pd.DataFrame(
            {
                "activity_id": [1, 2, 3, 4, 5, 6],
                "molecule_chembl_id": ["CHEMBL1"] * 6,
                "target_chembl_id": ["TARGET1"] * 6,
                "standard_smiles": ["CCO"] * 6,  # Same compound
                "canonical_smiles": ["CCO"] * 6,
                "pchembl_value": [6.0, 6.0, 6.0, 5.0, 7.0, 7.0],
                "standard_relation": ["<", "<", ">", "<", "=", "="],
                "mutation": ["WT"] * 6,
                "assay_chembl_id": ["ASSAY1"] * 6,
                "standard_type": ["IC50"] * 6,
                "assay_description": ["Test assay"] * 6,
                "assay_type": ["B"] * 6,
                "confidence_score": [9] * 6,
                "target_organism": ["Homo sapiens"] * 6,
                "document_chembl_id": ["DOC1"] * 6,
                "assay_tissue": [""] * 6,
                "assay_cell_type": [""] * 6,
                "relationship_type": ["D"] * 6,
                "max_phase": [1] * 6,
                "oral": [False] * 6,
                "prodrug": [False] * 6,
                "withdrawn_flag": [False] * 6,
                "doc_type": ["PUBLICATION"] * 6,
                "doi": ["10.1234/test"] * 6,
                "journal": ["Test Journal"] * 6,
                "year": [2020] * 6,
                "chembl_release": [30] * 6,
                "data_dropping_comment": [""] * 6,
                "data_processing_comment": [""] * 6,
            }
        )

        aggr_df = aggregate_data(
            test_data,
            chirality=False,
            metadata_cols=[],
            extra_id_cols=[],
            output_path=None,
            compound_equality="connectivity",
        )

        # Expected behavior:
        # Row 1,2: pchembl=6.0, relation='<' → should aggregate (same value AND relation)
        # Row 3: pchembl=6.0, relation='>' → separate row (different relation)
        # Row 4: pchembl=5.0, relation='<' → separate row (different value)
        # Row 5,6: pchembl=7.0, relation='=' → should aggregate (both are discrete measurements)

        # Check that we get 4 rows in the output
        self.assertEqual(len(aggr_df), 4, f"Expected 4 rows, got {len(aggr_df)}")

        # Check that rows with '<' and pchembl=6.0 were aggregated
        censored_lt_6 = aggr_df[
            (aggr_df["standard_relation"].str.contains("<", regex=False, na=False))
            & (aggr_df["pchembl_value_mean"].round(1) == 6.0)
        ]
        self.assertEqual(len(censored_lt_6), 1, "Expected 1 aggregated row for '<6.0' measurements")
        self.assertEqual(
            censored_lt_6.iloc[0]["pchembl_value_counts"], 2, "Expected 2 measurements aggregated for '<6.0'"
        )
        # For censored measurements, mean/median should be the unique value (not calculated)
        self.assertEqual(censored_lt_6.iloc[0]["pchembl_value_mean"], 6.0, "Expected mean to be 6.0 (unique value)")
        self.assertEqual(censored_lt_6.iloc[0]["pchembl_value_median"], 6.0, "Expected median to be 6.0 (unique value)")
        self.assertTrue(pd.isna(censored_lt_6.iloc[0]["pchembl_value_std"]), "Expected std to be NaN for censored data")

        # Check that row with '>' and pchembl=6.0 is separate
        censored_gt_6 = aggr_df[
            (aggr_df["standard_relation"].str.contains(">", regex=False, na=False))
            & (aggr_df["pchembl_value_mean"].round(1) == 6.0)
        ]
        self.assertEqual(len(censored_gt_6), 1, "Expected 1 row for '>6.0' measurement")
        self.assertEqual(
            censored_gt_6.iloc[0]["pchembl_value_counts"], 1, "Expected 1 count for single measurement '>6.0'"
        )

        # Check that row with '<' and pchembl=5.0 is separate
        censored_lt_5 = aggr_df[
            (aggr_df["standard_relation"].str.contains("<", regex=False, na=False))
            & (aggr_df["pchembl_value_mean"].round(1) == 5.0)
        ]
        self.assertEqual(len(censored_lt_5), 1, "Expected 1 row for '<5.0' measurement")
        self.assertEqual(
            censored_lt_5.iloc[0]["pchembl_value_counts"], 1, "Expected 1 count for single measurement '<5.0'"
        )

        # Check that discrete measurements with '=' were aggregated normally
        discrete_7 = aggr_df[
            (aggr_df["standard_relation"] == "=") & (aggr_df["pchembl_value_mean"].round(1) == 7.0)
        ]
        self.assertEqual(len(discrete_7), 1, "Expected 1 aggregated row for '=7.0' measurements")
        self.assertEqual(
            discrete_7.iloc[0]["pchembl_value_counts"], 2, "Expected 2 measurements aggregated for '=7.0'"
        )


if __name__ == "__main__":
    unittest.main()
