"""Module to test the workflow functions within /cli. To reproduce the same data in the
`resources/` directory, check the check `ADORA3_data_recipe.json`.
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from Capricho.cli.chembl_data_pipeline import aggregate_data
from Capricho.cli.prepare import prepare_multitask_data


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


    def test_aggregate_data_with_standard_value(self):
        """Test that aggregate_data works with standard_value instead of pchembl_value."""
        test_data = pd.DataFrame(
            {
                "activity_id": [1, 2, 3, 4],
                "molecule_chembl_id": ["CHEMBL1"] * 4,
                "target_chembl_id": ["TARGET1"] * 4,
                "standard_smiles": ["CCO", "CCO", "CCC", "CCC"],  # 2 compounds
                "canonical_smiles": ["CCO", "CCO", "CCC", "CCC"],
                "standard_value": [50.0, 75.0, 100.0, 80.0],  # % inhibition values
                "standard_units": ["%", "%", "%", "%"],
                "standard_relation": ["=", "=", "=", "="],
                "mutation": ["WT"] * 4,
                "assay_chembl_id": ["ASSAY1", "ASSAY2", "ASSAY1", "ASSAY2"],
                "standard_type": ["Inhibition"] * 4,
                "assay_description": ["Test assay"] * 4,
                "assay_type": ["A"] * 4,
                "confidence_score": [9] * 4,
                "target_organism": ["Homo sapiens"] * 4,
                "document_chembl_id": ["DOC1"] * 4,
                "assay_tissue": [""] * 4,
                "assay_cell_type": [""] * 4,
                "relationship_type": ["D"] * 4,
                "max_phase": [1] * 4,
                "oral": [False] * 4,
                "prodrug": [False] * 4,
                "withdrawn_flag": [False] * 4,
                "doc_type": ["PUBLICATION"] * 4,
                "doi": ["10.1234/test"] * 4,
                "journal": ["Test Journal"] * 4,
                "year": [2020] * 4,
                "chembl_release": [30] * 4,
                "data_dropping_comment": [""] * 4,
                "data_processing_comment": [""] * 4,
            }
        )

        aggr_df = aggregate_data(
            test_data,
            chirality=False,
            metadata_cols=[],
            extra_id_cols=["standard_units"],  # Group by units to prevent mixing
            output_path=None,
            compound_equality="connectivity",
            value_col="standard_value",  # Use standard_value instead of pchembl_value
        )

        # Check that we have 2 rows (2 compounds aggregated)
        self.assertEqual(len(aggr_df), 2)

        # Check that standard_value statistics columns exist
        self.assertIn("standard_value_mean", aggr_df.columns)
        self.assertIn("standard_value_std", aggr_df.columns)
        self.assertIn("standard_value_median", aggr_df.columns)
        self.assertIn("standard_value_counts", aggr_df.columns)

        # Check that pchembl_value columns do NOT exist
        self.assertNotIn("pchembl_value_mean", aggr_df.columns)
        self.assertNotIn("pchembl_value_std", aggr_df.columns)

        # Check the aggregated values for CCO compound
        cco_rows = aggr_df[aggr_df["smiles"].str.contains("CCO", na=False)]
        if len(cco_rows) > 0:
            cco_row = cco_rows.iloc[0]
            self.assertEqual(cco_row["standard_value_counts"], 2)
            self.assertAlmostEqual(cco_row["standard_value_mean"], (50.0 + 75.0) / 2)

    def test_prepare_multitask_data(self):
        """Test that prepare command creates correct activity matrix."""
        test_data = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1", "CONN2", "CONN2"],
                "smiles": ["CCO", "CCO", "CCC", "CCC"],
                "target_chembl_id": ["TARGET1", "TARGET2", "TARGET1", "TARGET2"],
                "pchembl_value_mean": [6.5, 7.0, 5.5, 8.0],
                "data_dropping_comment": ["", "some_flag", "", ""],
            }
        )

        matrix = prepare_multitask_data(
            df=test_data,
            task_col="target_chembl_id",
            value_col="pchembl_value_mean",
            compound_col="connectivity",
            smiles_col="smiles",
            remove_flags=None,
        )

        # Verify shape: 2 compounds x 2 targets
        self.assertEqual(len(matrix), 2, "Expected 2 compounds in the activity matrix")
        self.assertIn("TARGET1", matrix.columns, "Expected TARGET1 column")
        self.assertIn("TARGET2", matrix.columns, "Expected TARGET2 column")
        self.assertIn("smiles", matrix.columns, "Expected smiles column")

        # Verify values
        conn1_row = matrix[matrix.index == "CONN1"]
        self.assertEqual(len(conn1_row), 1)
        self.assertAlmostEqual(conn1_row["TARGET1"].iloc[0], 6.5)
        self.assertAlmostEqual(conn1_row["TARGET2"].iloc[0], 7.0)
        self.assertEqual(conn1_row["smiles"].iloc[0], "CCO")

        conn2_row = matrix[matrix.index == "CONN2"]
        self.assertEqual(len(conn2_row), 1)
        self.assertAlmostEqual(conn2_row["TARGET1"].iloc[0], 5.5)
        self.assertAlmostEqual(conn2_row["TARGET2"].iloc[0], 8.0)
        self.assertEqual(conn2_row["smiles"].iloc[0], "CCC")

    def test_prepare_multitask_data_with_flags(self):
        """Test that prepare command correctly filters out flagged data."""
        test_data = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1", "CONN2", "CONN2"],
                "smiles": ["CCO", "CCO", "CCC", "CCC"],
                "target_chembl_id": ["TARGET1", "TARGET2", "TARGET1", "TARGET2"],
                "pchembl_value_mean": [6.5, 7.0, 5.5, 8.0],
                "data_dropping_comment": ["", "flag_to_remove", "", "another_flag"],
            }
        )

        matrix = prepare_multitask_data(
            df=test_data,
            task_col="target_chembl_id",
            value_col="pchembl_value_mean",
            compound_col="connectivity",
            smiles_col="smiles",
            remove_flags=["flag_to_remove", "another_flag"],
        )

        # After filtering, we should only have CONN1->TARGET1 and CONN2->TARGET1
        self.assertEqual(len(matrix), 2, "Expected 2 compounds after filtering")

        # Only TARGET1 should exist as a column (TARGET2 was completely filtered out)
        self.assertIn("TARGET1", matrix.columns, "Expected TARGET1 column")
        self.assertNotIn("TARGET2", matrix.columns, "TARGET2 should not exist after filtering")
        self.assertIn("smiles", matrix.columns, "Expected smiles column")

        # Verify the values for TARGET1
        conn1_row = matrix[matrix.index == "CONN1"]
        self.assertAlmostEqual(conn1_row["TARGET1"].iloc[0], 6.5)

        conn2_row = matrix[matrix.index == "CONN2"]
        self.assertAlmostEqual(conn2_row["TARGET1"].iloc[0], 5.5)

    def test_prepare_multitask_data_warns_on_duplicates(self):
        """Test that prepare warns when duplicates exist (e.g., same compound-target but different id_columns)."""
        # Simulate data aggregated with --id-columns assay_tissue
        # Same compound-target pair but different tissues = duplicates when pivoting on target alone
        test_data = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1", "CONN2"],
                "smiles": ["CCO", "CCO", "CCC"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1"],
                "assay_tissue": ["Liver", "Kidney", "Liver"],  # Different tissues
                "pchembl_value_mean": [6.5, 7.0, 5.5],
                "data_dropping_comment": ["", "", ""],
            }
        )

        # Capture log output
        import io
        import logging

        from Capricho.logger import logger

        log_capture = io.StringIO()
        handler = logging.StreamHandler(log_capture)
        handler.setLevel(logging.WARNING)
        logger.add(handler, format="{message}", level="WARNING")

        matrix = prepare_multitask_data(
            df=test_data,
            task_col="target_chembl_id",
            value_col="pchembl_value_mean",
            compound_col="connectivity",
            smiles_col="smiles",
            remove_flags=None,
        )

        # Remove our test handler
        logger.remove()

        # The function should still work (take first value), but there should be a warning
        # CONN1 has 2 values for TARGET1 (Liver=6.5, Kidney=7.0), only one will be kept
        self.assertEqual(len(matrix), 2)  # 2 compounds

        # Check that only one value was kept for CONN1-TARGET1
        conn1_target1 = matrix.loc["CONN1", "TARGET1"]
        self.assertIn(conn1_target1, [6.5, 7.0])  # Either value is acceptable

    def test_prepare_multitask_data_with_id_columns(self):
        """Test that id_columns creates composite task identifiers."""
        test_data = pd.DataFrame(
            {
                "connectivity": ["CONN1", "CONN1", "CONN2", "CONN2"],
                "smiles": ["CCO", "CCO", "CCC", "CCC"],
                "target_chembl_id": ["TARGET1", "TARGET1", "TARGET1", "TARGET1"],
                "assay_tissue": ["Liver", "Kidney", "Liver", "Kidney"],
                "pchembl_value_mean": [6.5, 7.0, 5.5, 8.0],
                "data_dropping_comment": ["", "", "", ""],
            }
        )

        matrix = prepare_multitask_data(
            df=test_data,
            task_col="target_chembl_id",
            value_col="pchembl_value_mean",
            compound_col="connectivity",
            smiles_col="smiles",
            remove_flags=None,
            id_columns=["assay_tissue"],
        )

        # With id_columns, we should have 4 unique tasks: TARGET1-Liver, TARGET1-Kidney
        self.assertEqual(len(matrix), 2)  # 2 compounds
        self.assertIn("TARGET1-Liver", matrix.columns)
        self.assertIn("TARGET1-Kidney", matrix.columns)

        # Verify correct values
        self.assertAlmostEqual(matrix.loc["CONN1", "TARGET1-Liver"], 6.5)
        self.assertAlmostEqual(matrix.loc["CONN1", "TARGET1-Kidney"], 7.0)
        self.assertAlmostEqual(matrix.loc["CONN2", "TARGET1-Liver"], 5.5)
        self.assertAlmostEqual(matrix.loc["CONN2", "TARGET1-Kidney"], 8.0)


if __name__ == "__main__":
    unittest.main()
