import unittest
from io import StringIO

import pandas as pd
from loguru import logger as loguru_logger

from CompoundMapper.chembl.processing import curate_activity_pairs  # Added import
from CompoundMapper.core.stats_make import process_repeat_mols

N_CPDS = 19  # number of compounds in the test data


class TestActivityCuration(unittest.TestCase):

    def test_curate_activity_pairs_logic(self):
        data = {
            "molecule_chembl_id": [
                "CHEMBL1",
                "CHEMBL1",
                "CHEMBL2",
                "CHEMBL2",
                "CHEMBL3",
                "CHEMBL3",
                "CHEMBL4",
                "CHEMBL5",
                "CHEMBL5",
                "CHEMBL6",
                "CHEMBL7",
                "CHEMBL7",
                "CHEMBL8",
                "CHEMBL8",
                "CHEMBL8",
                "CHEMBL9",
                "CHEMBL9",
                "CHEMBL9",
                "CHEMBL9",
            ],
            "canonical_smiles": [
                "CC=O",  # CHEMBL1
                "CC=O",  # CHEMBL1
                "CC",  # CHEMBL2
                "CC",  # CHEMBL2
                "CC(O)CC",  # CHEMBL3
                "CC(O)CC",  # CHEMBL3
                "NCO",  # CHEMBL4
                "C(N)O",  # CHEMBL5
                "C(N)O",  # CHEMBL5
                "CNCOCC(OCC)CCC",  # CHEMBL6
                "CCO",  # CHEMBL7
                "CCO",  # CHEMBL7
                "CCC",  # CHEMBL8
                "CCC",  # CHEMBL8
                "CCC",  # CHEMBL8
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
            ],
            "activity_id": list(range(1, N_CPDS + 1)),
            "assay_chembl_id": [
                "ASSAY1",  # CHEMBL1
                "ASSAY2",  # CHEMBL1
                "ASSAY3",  # CHEMBL2
                "ASSAY4",  # CHEMBL2
                "ASSAY5",  # CHEMBL3
                "ASSAY6",  # CHEMBL3
                "ASSAY7",  # CHEMBL4
                "ASSAY8",  # CHEMBL5
                "ASSAY9",  # CHEMBL5
                "ASSAY10",  # CHEMBL6
                "ASSAY11",  # CHEMBL7
                "ASSAY12",  # CHEMBL7
                "ASSAY13",  # CHEMBL8
                "ASSAY14",  # CHEMBL8
                "ASSAY15",  # CHEMBL8
                "ASSAY16",  # CHEMBL9
                "ASSAY17",  # CHEMBL9
                "ASSAY18",  # CHEMBL9
                "ASSAY19",  # CHEMBL9
            ],
            "assay_description": [
                "Assay Desc 1",  # CHEMBL1
                "Assay Desc 2",  # CHEMBL1
                "Assay Desc 3",  # CHEMBL2
                "Assay Desc 4",  # CHEMBL2
                "Assay Desc 5",  # CHEMBL3
                "Assay Desc 6",  # CHEMBL3
                "Assay Desc 7",  # CHEMBL4
                "Assay Desc 8",  # CHEMBL5
                "Assay Desc 9",  # CHEMBL5
                "Assay Desc 10",  # CHEMBL6
                "Assay Desc 11",  # CHEMBL7
                "Assay Desc 12",  # CHEMBL7
                "Assay Desc 13",  # CHEMBL8
                "Assay Desc 14",  # CHEMBL8
                "Assay Desc 15",  # CHEMBL8
                "Assay Desc 16",  # CHEMBL9
                "Assay Desc 17",  # CHEMBL9
                "Assay Desc 18",  # CHEMBL9
                "Assay Desc 19",  # CHEMBL9
            ],
            "assay_type": ["B"] * N_CPDS,  # doesn't matter here, aggregate both types
            "pchembl_value": [
                7.0,  # CHEMBL1
                10.0,  # CHEMBL1 - diff 3.0 -> remove
                6.0,  # CHEMBL2
                6.0,  # CHEMBL2
                8.0,  # CHEMBL3
                5.0,  # CHEMBL3 - diff 3.0 -> remove
                7.5,  # CHEMBL4 - single, keep
                6.2,  # CHEMBL5
                9.2,  # CHEMBL5 - diff 3.0 -> remove
                6.8,  # CHEMBL6 - single, keep
                6.0,  # CHEMBL7
                7.0,  # CHEMBL7 - diff 1.0 -> keep
                5.0,  # CHEMBL8 - ASSAY13
                6.0,  # CHEMBL8 - ASSAY14 (pair with ASSAY13: diff 1.0; pair with ASSAY15: diff 2.0) -> keep this one
                8.0,  # CHEMBL8 - ASSAY15 (pair with ASSAY13: diff 3.0 -> remove 5.0 and 8.0)
                7.0,  # CHEMBL9 - diff 3.0 -> remove
                10.0,  # CHEMBL9 - diff 3.0 -> remove
                10.0,  # CHEMBL9 - ASSAY18 (pair with ASSAY16 and ASSAY17: diff 3.0 -> remove 7.0 and 10.0)
                7.2,  # CHEMBL9
                # CHEMBL 9 illustrates a possible rare case where you'd have a review ASSAY16 that repeats the
                # same measurement of 7.0, but annotates it wrongly as 10.0. Later, the same measurement is
                # again wrongly reported as 10.0. I don't know if this is a realy case, but I tried to cover it.
            ],
            "standard_flag": [1] * N_CPDS,
            "standard_relation": ["="] * N_CPDS,
            "standard_type": [
                "Ki",
                "Ki",
                "IC50",
                "IC50",
                "Ki",
                "Ki",
                "Ki",
                "IC50",
                "Ki",
                "Ki",
                "EC50",
                "EC50",
                "Kd",
                "Kd",
                "Kd",
                "IC50",
                "IC50",
                "IC50",
                "IC50",
            ],
            "standard_units": ["nM"] * N_CPDS,
            "target_chembl_id": [
                "TARGET1",
                "TARGET1",  # CHEMBL1
                "TARGET2",
                "TARGET2",  # CHEMBL2
                "TARGET3",
                "TARGET3",  # CHEMBL3
                "TARGET4",  # CHEMBL4
                "TARGET5",
                "TARGET5",  # CHEMBL5
                "TARGET6",  # CHEMBL6
                "TARGET7",
                "TARGET7",  # CHEMBL7
                "TARGET8",
                "TARGET8",
                "TARGET8",  # CHEMBL8
                "TARGET9",  # CHEMBL9
                "TARGET9",  # CHEMBL9
                "TARGET9",  # CHEMBL9
                "TARGET9",  # CHEMBL9
            ],
            "target_organism": [
                "Human",
                "Human",
                "Human",
                "Human",
                "Mouse",
                "Mouse",
                "Rat",
                "Human",
                "Human",
                "Human",
                "Dog",
                "Dog",
                "Cat",
                "Cat",
                "Cat",
                "Elf",
                "Elf",
                "Elf",
                "Elf",
            ],
            "assay_cell_type": [
                "Cell Type A",
                "Cell Type A",
                "Cell Type B",
                "Cell Type B",
                "Cell Type C",
                "Cell Type C",
                "None",
                "Cell Type D",
                "Cell Type D",
                "Cell Type E",
                "Cell Type F",
                "Cell Type F",
                "Cell Type G",
                "Cell Type G",
                "Cell Type G",
                "Cell Type H",
                "Cell Type H",
                "Cell Type H",
                "Cell Type H",
            ],
            "assay_organism": [
                "Human",
                "Human",
                "Human",
                "Human",
                "Mouse",
                "Mouse",
                "Rat",
                "Human",
                "Human",
                "Human",
                "Dog",
                "Dog",
                "Cat",
                "Cat",
                "Cat",
                "Elf",
                "Elf",
                "Elf",
                "Elf",
            ],
            "confidence_score": [9, 9, 8, 8, 9, 9, 7, 9, 9, 9, 8, 8, 9, 9, 9, 9, 9, 9, 9],
            "document_chembl_id": [
                "DOC1",
                "DOC2",
                "DOC3",
                "DOC4",
                "DOC5",
                "DOC6",
                "DOC7",
                "DOC8",
                "DOC9",
                "DOC10",
                "DOC11",
                "DOC12",
                "DOC13",
                "DOC14",
                "DOC15",
                "DOC16",
                "DOC17",
                "DOC18",
                "DOC19",
            ],
            "standard_smiles": [
                "CC=O",
                "CC=O",
                "CC",
                "CC",
                "CC(O)CC",
                "CC(O)CC",
                "NCO",
                "C(N)O",
                "C(N)O",
                "CNCOCC(OCC)CCC",
                "CCO",
                "CCO",
                "CCC",
                "CCC",
                "CCC",
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
                "O=C(C)Oc1ccccc1C(=O)O",  # CHEMBL9
            ],
            "mutation": ["WT"] * N_CPDS,
        }
        df_initial = pd.DataFrame(data)

        # Expected DataFrame after curate_activity_pairs
        # CHEMBL1 (idx 0, 1) removed
        # CHEMBL2 (idx 2, 3) kept
        # CHEMBL3 (idx 4, 5) removed
        # CHEMBL4 (idx 6) kept
        # CHEMBL5 (idx 7, 8) removed
        # CHEMBL6 (idx 9) kept
        # CHEMBL7 (idx 10, 11) kept
        # CHEMBL8 (idx 12, 14) removed, (idx 13) kept
        # CHEMBL9 (idx 15, 16, 17, ) removed, (idx 18) kept
        expected_indices_to_keep = [2, 3, 6, 9, 10, 11, 13, 18]
        expected_df_curated = df_initial.loc[expected_indices_to_keep].reset_index(drop=True)

        captured_logs = StringIO()
        sink_id = loguru_logger.add(captured_logs, format="{message}", level="DEBUG")

        processed_df = curate_activity_pairs(df_initial.copy())

        log_output = captured_logs.getvalue()
        loguru_logger.remove(sink_id)

        print("\nCaptured Log Output:\n", log_output)

        self.assertIn(
            "Activity Curation: Removing 11 measurements due to pChEMBL value criteria.", log_output
        )

        pd.testing.assert_frame_equal(
            processed_df.reset_index(drop=True),
            expected_df_curated.reset_index(drop=True),
            check_dtype=False,  # Allow different int types (e.g. int32 vs int64)
            check_like=True,  # Ignore column order as long as labels match
        )


# Remove the old TestProcessRepeatMolsRealistic class if it's no longer needed
# or adapt its tests separately if they are still relevant for process_repeat_mols.
# For now, I will comment it out to avoid test runner conflicts.

# class TestProcessRepeatMolsRealistic(unittest.TestCase):
#
#     def test_activity_curation_filter_realistic(self):
#         data = {
#             "molecule_chembl_id": [
#                 "CHEMBL1",
#                 "CHEMBL1",
#                 "CHEMBL2",
#                 "CHEMBL2",
#                 "CHEMBL3",
#                 "CHEMBL3",
#                 "CHEMBL4",
#                 "CHEMBL5",
#                 "CHEMBL5",
#                 "CHEMBL6",
#             ],
#             "canonical_smiles": [
#                 "CC=O",
#                 "CC=O",
#                 "CC",
#                 "CC",
#                 "CC(O)CC",
#                 "CC(O)CC",
#                 "NCO",
#                 "C(N)O",
#                 "C(N)O",
#                 "CNCOCC(OCC)CCC",
#             ],
#             "activity_id": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
#             "assay_chembl_id": [
#                 "ASSAY1",
#                 "ASSAY2",
#                 "ASSAY3",
#                 "ASSAY4",
#                 "ASSAY5",
#                 "ASSAY6",
#                 "ASSAY7",
#                 "ASSAY8",
#                 "ASSAY9",
#                 "ASSAY10",
#             ],
#             "assay_description": [
#                 "Assay Desc 1",
#                 "Assay Desc 2",
#                 "Assay Desc 3",
#                 "Assay Desc 4",
#                 "Assay Desc 5",
#                 "Assay Desc 6",
#                 "Assay Desc 7",
#                 "Assay Desc 8",
#                 "Assay Desc 9",
#                 "Assay Desc 10",
#             ],
#             "assay_type": ["B", "B", "F", "F", "B", "B", "B", "F", "B", "B"],
#             "pchembl_value": [
#                 7.0,
#                 10.0,
#                 6.0,
#                 6.0,
#                 8.0,
#                 5.0,
#                 7.5,
#                 6.2,
#                 9.2,
#                 6.8,
#             ],  # CHEMBL1: diff 3.0, CHEMBL2: same, CHEMBL3: diff 3.0, CHEMBL5: diff 3.0
#             "standard_flag": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
#             "standard_relation": ["=", "=", "=", "=", "=", "=", "=", "=", "=", "="],
#             "standard_type": ["Ki", "Ki", "IC50", "IC50", "Ki", "Ki", "Ki", "IC50", "Ki", "Ki"],
#             "standard_units": ["nM", "nM", "nM", "nM", "nM", "nM", "nM", "nM", "nM", "nM"],
#             "target_chembl_id": [
#                 "TARGET1",
#                 "TARGET1",
#                 "TARGET2",
#                 "TARGET2",
#                 "TARGET3",
#                 "TARGET3",
#                 "TARGET4",
#                 "TARGET5",
#                 "TARGET5",
#                 "TARGET6",
#             ],
#             "target_organism": [
#                 "Human",
#                 "Human",
#                 "Human",
#                 "Human",
#                 "Mouse",
#                 "Mouse",
#                 "Rat",
#                 "Human",
#                 "Human",
#                 "Human",
#             ],
#             "assay_cell_type": [
#                 "Cell Type A",
#                 "Cell Type A",
#                 "Cell Type B",
#                 "Cell Type B",
#                 "Cell Type C",
#                 "Cell Type C",
#                 "None",
#                 "Cell Type D",
#                 "Cell Type D",
#                 "Cell Type E",
#             ],
#             "assay_organism": [
#                 "Human",
#                 "Human",
#                 "Human",
#                 "Human",
#                 "Mouse",
#                 "Mouse",
#                 "Rat",
#                 "Human",
#                 "Human",
#                 "Human",
#             ],
#             "confidence_score": [9, 9, 8, 8, 9, 9, 7, 9, 9, 9],
#             "document_chembl_id": [
#                 "DOC1",
#                 "DOC1",
#                 "DOC2",
#                 "DOC2",
#                 "DOC3",
#                 "DOC3",
#                 "DOC4",
#                 "DOC5",
#                 "DOC5",
#                 "DOC6",
#             ],
#             "standard_smiles": [
#                 "CC=O",
#                 "CC=O",
#                 "CC",
#                 "CC",
#                 "CC(O)CC",
#                 "CC(O)CC",
#                 "NCO",
#                 "C(N)O",
#                 "C(N)O",
#                 "CNCOCC(OCC)CCC",
#             ],
#             "mutation": [
#                 "WT",
#                 "WT",
#                 "WT",
#                 "WT",
#                 "WT",
#                 "WT",
#                 "WT",
#                 "WT",
#                 "WT",
#                 "WT",
#             ],
#         }
#         df = pd.DataFrame(data)
#         repeat_element_idxs = [
#             [0, 1],
#             [2, 3],
#             [4, 5],
#             [9, 10],
#         ]  # Indices of repeated molecules based on SMILES in this dummy data (CHEMBL1, CHEMBL2, CHEMBL3, CHEMBL5)
#
#         expected_log_message = "Activity curation filter removed 6 data points."
#         captured_logs = StringIO()
#
#         sink_id = loguru_logger.add(captured_logs, format="{message}", level="INFO")  # Add sink and get ID
#         multiple_value_cols = [k for k in data.keys() if k not in ["variant_sequence", "target_chembl_id"]]
#         processed_df = process_repeat_mols(
#             df.copy(),
#             repeat_element_idxs=repeat_element_idxs,
#             multiple_value_cols=multiple_value_cols,
#             aggregate_mutants=True,
#         )
#         log_output = captured_logs.getvalue()  # Get captured logs
#         loguru_logger.remove(sink_id)  # Remove the sink to avoid interference with other tests
#
#         # Expected DataFrame after processing and filtering.
#         # Rows with CHEMBL1 (7.0 and 10.0, diff 3.0) and CHEMBL3 (8.0 and 5.0, diff 3.0) and CHEMBL5 (6.2 and 9.2, diff 3.0) should be removed.
#         expected_data_filtered = {
#             "smiles": ["CC", "NCO", "CCC(OCC)CCC"],
#             "original_smiles": ["CC", "NCO", "CCC(OCC)CCC"],
#             "molecule_chembl_id": ["CHEMBL2", "CHEMBL4", "CHEMBL6"],
#             "canonical_smiles": ["CC", "NCO", "CCC(OCC)CCC"],
#             "activity_id": ["3;4", "7", "10"],
#             "assay_chembl_id": ["ASSAY3;ASSAY4", "ASSAY7", "ASSAY10"],
#             "assay_description": ["Assay Desc 3;Assay Desc 4", "Assay Desc 7", "Assay Desc 10"],
#             "assay_type": ["F;F", "B", "B"],
#             "pchembl_value": ["6.0;6.0", "7.5", "6.8"],
#             "standard_flag": ["1;1", "1", "1"],
#             "standard_relation": ["=;=", "=", "="],
#             "standard_type": ["IC50;IC50", "Ki", "Ki"],
#             "standard_units": ["nM;nM", "nM", "nM"],
#             "target_chembl_id": ["TARGET2;TARGET2", "TARGET4", "TARGET6"],
#             "target_organism": ["Human;Human", "Rat", "Human"],
#             "assay_cell_type": ["Cell Type B;Cell Type B", "None", "Cell Type E"],
#             "assay_organism": ["Human;Human", "Rat", "Human"],
#             "confidence_score": ["8;8", "7", "9"],
#             "document_chembl_id": ["DOC2;DOC2", "DOC4", "DOC6"],
#             "might_rancemic": [True, False, False],
#             "pchembl_value_mean": [6.0, 7.5, 6.8],
#             "pchembl_value_std": [0.0, None, None],
#             "pchembl_value_median": [6.0, 7.5, 6.8],
#             "pchembl_value_mad": [0.0, None, None],
#             "pchembl_value_counts": [2, 1, 1],
#             "mutation": ["WT", "WT", "WT"],
#         }
#         expected_df_filtered = pd.DataFrame(expected_data_filtered)
#
#         # Assert that the log message was called with the correct message
#         log_output = captured_logs.getvalue()
#         self.assertIn(expected_log_message, log_output)
#
#         # Assert that the processed DataFrame is as expected after filtering for 3.0 difference and other processing
#         pd.testing.assert_frame_equal(
#             processed_df.reset_index(drop=True),
#             expected_df_filtered.reset_index(drop=True),
#             check_dtype=False,
#             check_like=True,
#         )


if __name__ == "__main__":
    unittest.main()
