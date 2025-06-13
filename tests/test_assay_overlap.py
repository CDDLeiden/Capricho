import unittest
from io import StringIO

import pandas as pd

from CompoundMapper.chembl.data_flag_functions import flag_insufficient_assay_overlap
from CompoundMapper.core.default_fields import (
    ACTIVITY_ID,
    ASSAY_ID,
    DATA_DROPPING_COMMENT,
    MOLECULE_ID,
    TARGET_ID,
)
from CompoundMapper.logger import logger


class TestFlagInsufficientAssayOverlap(unittest.TestCase):

    def setUp(self):
        self.activity_col = ACTIVITY_ID  # Unique identifier for each activity
        self.molecule_col = MOLECULE_ID
        self.assay_col = ASSAY_ID
        self.target_col = TARGET_ID
        self.comment_col = DATA_DROPPING_COMMENT

        test_scenario_data = [
            # Target T1:l
            # A1: M1, M2, M3
            # A2: M1, M4
            # Overlap(A1,A2) for T1 = 1 (M1)
            ("T1", "A1", [f"M{i}" for i in range(1, 4)]),
            ("T1", "A2", ["M1", "M4"]),
            # Target T2:
            # A3: M10, M11, M12
            # A4: M10, M11, M13
            # A5: M10, M14
            # Overlap(A3,A4) for T2 = 2 (M10, M11)
            # Overlap(A3,A5) for T2 = 1 (M10)
            # Overlap(A4,A5) for T2 = 1 (M10)
            ("T2", "A3", [f"M{i}" for i in range(10, 13)]),
            ("T2", "A4", ["M10", "M11", "M13"]),
            ("T2", "A5", ["M10", "M14"]),
            # Target T3:
            # A6: M20, M21, M22, M23, M24
            # A7: M20, M21, M22, M25, M26
            # Overlap(A6,A7) for T3 = 3 (M20, M21, M22)
            ("T3", "A6", [f"M{i}" for i in range(20, 25)]),
            ("T3", "A7", [f"M{i}" for i in range(20, 23)] + [f"M{i}" for i in range(25, 27)]),
            # Target T4: Single assay A8
            ("T4", "A8", [f"M{i}" for i in range(30, 32)]),
        ]

        data_for_df = []
        activity_id_counter = 1
        for target, assay, molecules in test_scenario_data:
            for molecule in molecules:
                data_for_df.append(
                    {
                        self.target_col: target,
                        self.assay_col: assay,
                        self.molecule_col: molecule,
                        # unique identifier / activity like chembl
                        self.activity_col: activity_id_counter,
                        self.comment_col: None,
                    }
                )
                activity_id_counter += 1

        self.df_initial = pd.DataFrame(data_for_df)

    def tearDown(self):
        # General cleanup if any sinks were missed, though specific tests should handle their own.
        # This might not be strictly necessary if tests are isolated.
        pass

    def _run_test_and_asserts(
        self,
        df_input,
        min_overlap,
        expected_flagged_assays_global,
        expected_total_flagged_activities,
        expected_total_flagged_unique_assays,
        check_no_overlap_log=False,
    ):
        logs_io = StringIO()
        sink_id = logger.add(logs_io, format="{message}", level="DEBUG", enqueue=True)

        try:
            df_processed = flag_insufficient_assay_overlap(
                df_input.copy(),
                min_overlap=min_overlap,
                molecule_col=self.molecule_col,
                assay_col=self.assay_col,
                target_col=self.target_col,
                comment_col=self.comment_col,
            )

            df_expected = df_input.copy()
            if expected_flagged_assays_global:
                expected_comment_text = f"Insufficient assay overlap (min_overlap={min_overlap})"
                mask = df_expected[self.assay_col].isin(list(expected_flagged_assays_global))
                df_expected.loc[mask, self.comment_col] = expected_comment_text

            logger.remove(sink_id)  # Remove sink before getting value to ensure logs are flushed
            sink_id = None  # Mark sink_id as removed
            log_output = logs_io.getvalue()
            print(f"\n--- Logs for min_overlap={min_overlap} ---\n{log_output}\n--- End Logs ---")

            if expected_total_flagged_activities > 0:
                self.assertIn(
                    f"Flagging {expected_total_flagged_activities} activities from {expected_total_flagged_unique_assays} unique assays "
                    f"that were part of pairs not meeting the minimum overlap of {min_overlap} compounds.",
                    log_output,
                )
            elif check_no_overlap_log:
                self.assertIn(
                    f"No assay pairs found below the minimum overlap of {min_overlap} compounds.", log_output
                )

            pd.testing.assert_frame_equal(
                df_processed.sort_values(by=[self.activity_col]).reset_index(drop=True),
                df_expected.sort_values(by=[self.activity_col]).reset_index(drop=True),
                check_dtype=False,
            )
        finally:
            if sink_id is not None and hasattr(logger._core, "handlers") and sink_id in logger._core.handlers:
                logger.remove(sink_id)
            logs_io.close()

    def test_min_overlap_flags_correctly(self):
        # Case 1: min_overlap = 2
        # T1: A1-A2 overlap 1. Flag A1, A2.
        # T2: A3-A4 overlap 2 (OK). A3-A5 overlap 1 (Flag A3,A5). A4-A5 overlap 1 (Flag A4,A5).
        # Globally flagged assays: A1, A2, A3, A4, A5
        # Activities: A1(3), A2(2), A3(3), A4(3), A5(2) = 13 activities from 5 unique assays
        self._run_test_and_asserts(self.df_initial, 2, {"A1", "A2", "A3", "A4", "A5"}, 13, 5)

    def test_min_overlap_all_pass(self):
        # Case 2: min_overlap = 1 (all pairs have at least 1 overlap)
        # No assays should be flagged.
        self._run_test_and_asserts(self.df_initial, 1, set(), 0, 0, check_no_overlap_log=True)

    def test_min_overlap_none_pass(self):
        # Case 3: min_overlap = 10 (no pairs will meet this)
        # T1: A1,A2 flagged.
        # T2: A3,A4,A5 flagged.
        # T3: A6,A7 flagged.
        # Globally flagged assays: A1, A2, A3, A4, A5, A6, A7
        # Activities: A1(3),A2(2) + A3(3),A4(3),A5(2) + A6(5),A7(5) = 5 + 8 + 10 = 23 activities from 7 unique assays
        self._run_test_and_asserts(self.df_initial, 10, {"A1", "A2", "A3", "A4", "A5", "A6", "A7"}, 23, 7)

    def test_min_overlap_is_zero(self):
        logs_io = StringIO()
        sink_id = logger.add(logs_io, format="{message}", level="DEBUG", enqueue=True)
        try:
            df_processed = flag_insufficient_assay_overlap(
                self.df_initial.copy(),
                min_overlap=0,
                molecule_col=self.molecule_col,
                assay_col=self.assay_col,
                target_col=self.target_col,
                comment_col=self.comment_col,
            )
            pd.testing.assert_frame_equal(df_processed, self.df_initial)
            logger.remove(sink_id)
            log_output = logs_io.getvalue()
            self.assertIn("min_overlap set to 0. Skipping insufficient assay overlap filtering.", log_output)
        finally:
            if hasattr(logger._core, "handlers") and sink_id in logger._core.handlers:
                logger.remove(sink_id)
            logs_io.close()

    def test_empty_dataframe(self):
        df_empty = pd.DataFrame(
            columns=[self.target_col, self.assay_col, self.molecule_col, self.comment_col, self.activity_col]
        )
        logs_io = StringIO()
        sink_id = logger.add(logs_io, format="{message}", level="DEBUG", enqueue=True)
        try:
            df_processed = flag_insufficient_assay_overlap(
                df_empty.copy(),
                min_overlap=2,
                molecule_col=self.molecule_col,
                assay_col=self.assay_col,
                target_col=self.target_col,
                comment_col=self.comment_col,
            )
            self.assertTrue(df_processed.empty)
            logger.remove(sink_id)
            log_output = logs_io.getvalue()
            self.assertIn("DataFrame is empty. Skipping insufficient assay overlap filtering.", log_output)
        finally:
            if hasattr(logger._core, "handlers") and sink_id in logger._core.handlers:
                logger.remove(sink_id)
            logs_io.close()

    def test_missing_columns(self):
        df_missing_col = self.df_initial.drop(columns=[self.target_col])
        logs_io = StringIO()
        sink_id = logger.add(logs_io, format="{message}", level="DEBUG", enqueue=True)
        try:
            df_processed = flag_insufficient_assay_overlap(
                df_missing_col.copy(),
                min_overlap=2,
                molecule_col=self.molecule_col,
                assay_col=self.assay_col,
                target_col=self.target_col,
                comment_col=self.comment_col,
            )
            logger.remove(sink_id)  # remove sink before getting value
            log_output = logs_io.getvalue()

            self.assertIn(
                f"One or more required columns ({self.molecule_col}, {self.assay_col}, {self.target_col}) "
                "not found in DataFrame. Skipping insufficient assay overlap filtering.",
                log_output,
            )
            # The function should return the original df (df_missing_col) in this case
            pd.testing.assert_frame_equal(df_processed, df_missing_col)
        finally:
            if hasattr(logger._core, "handlers") and sink_id in logger._core.handlers:
                logger.remove(sink_id)
            logs_io.close()


if __name__ == "__main__":
    unittest.main()
