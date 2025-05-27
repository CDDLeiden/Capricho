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
            "Activity Curation: Removing 11 measurements. potential unit annotation error "
            "- pChEMBL values differ by 3.0 for same molecule",
            log_output,
        )

        pd.testing.assert_frame_equal(
            processed_df.reset_index(drop=True),
            expected_df_curated.reset_index(drop=True),
            check_dtype=False,  # Allow different int types (e.g. int32 vs int64)
            check_like=True,  # Ignore column order as long as labels match
        )


if __name__ == "__main__":
    unittest.main()
