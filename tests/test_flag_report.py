"""Tests for the quality-flag accounting reported by the CLI and the case studies."""

import pandas as pd
import pytest

from Capricho.analysis import DroppingComment
from Capricho.flag_report import (
    ANY_FLAG_LABEL,
    ANY_SELECTED_LABEL,
    RETAINED_LABEL,
    TOTAL_LABEL,
    UNFLAGGED_LABEL,
    flags_in_comment,
    format_cross_assay_coverage,
    format_flag_summary,
    measurement_comments,
    summarize_cross_assay_coverage,
    summarize_flags,
)


def count_of(summary: pd.DataFrame, flag: str) -> int:
    """Absolute count reported for a flag."""
    return int(summary.loc[summary["flag"] == flag, "n"].iloc[0])


def pct_of(summary: pd.DataFrame, flag: str) -> float:
    """Percentage reported for a flag."""
    return float(summary.loc[summary["flag"] == flag, "pct"].iloc[0])


class TestMeasurementComments:
    def test_splits_pooled_measurements(self):
        assert measurement_comments("A|B|C") == ["A", "B", "C"]

    def test_single_comment_stays_whole(self):
        assert measurement_comments("Undefined Stereochemistry") == ["Undefined Stereochemistry"]

    def test_missing_values_yield_one_empty_comment(self):
        assert measurement_comments(None) == [""]
        assert measurement_comments(float("nan")) == [""]
        assert measurement_comments("nan") == [""]

    def test_empty_positions_are_kept_so_measurements_stay_aligned(self):
        assert measurement_comments("|Potential Duplicate|") == ["", "Potential Duplicate", ""]


class TestFlagsInComment:
    def test_splits_flags_joined_on_one_measurement(self):
        assert flags_in_comment("Potential Duplicate & Undefined Stereochemistry") == {
            "Potential Duplicate",
            "Undefined Stereochemistry",
        }

    def test_normalizes_run_specific_thresholds(self):
        assert flags_in_comment("Assay size < 5") == {"Assay size <"}
        assert flags_in_comment("Assay size < 20") == {"Assay size <"}

    def test_empty_comment_has_no_flags(self):
        assert flags_in_comment("") == set()
        assert flags_in_comment("nan") == set()


class TestSummarizeFlags:
    @pytest.fixture
    def flagged(self):
        """Four measurements: two share a flag, one carries two flags, one is clean."""
        return pd.DataFrame(
            {
                "data_dropping_comment": [
                    "Undefined Stereochemistry",
                    "Undefined Stereochemistry & Assay size < 5",
                    "Potential Duplicate",
                    "",
                ]
            }
        )

    def test_counts_each_flag_and_its_share(self, flagged):
        summary = summarize_flags(flagged)

        assert count_of(summary, "Undefined Stereochemistry") == 2
        assert pct_of(summary, "Undefined Stereochemistry") == 50.0
        assert count_of(summary, "Assay size <") == 1
        assert count_of(summary, "Potential Duplicate") == 1

    def test_union_counts_co_flagged_entries_once(self, flagged):
        summary = summarize_flags(flagged)

        # 2 + 1 + 1 = 4 per-flag hits, but only 3 distinct flagged measurements.
        assert summary.loc[summary["kind"] == "flag", "n"].sum() == 4
        assert count_of(summary, ANY_FLAG_LABEL) == 3
        assert pct_of(summary, ANY_FLAG_LABEL) == 75.0

    def test_totals_partition_the_dataset(self, flagged):
        summary = summarize_flags(flagged)

        assert count_of(summary, UNFLAGGED_LABEL) == 1
        assert count_of(summary, TOTAL_LABEL) == len(flagged)
        assert count_of(summary, ANY_FLAG_LABEL) + count_of(summary, UNFLAGGED_LABEL) == len(flagged)

    def test_complement_is_named_for_a_filter_only_when_flags_are_selected(self, flagged):
        """Without a selection the summary describes the data; with one it describes a filter."""
        assert UNFLAGGED_LABEL in set(summarize_flags(flagged)["flag"])
        assert RETAINED_LABEL in set(summarize_flags(flagged, flags=["Potential Duplicate"])["flag"])

    def test_selected_flags_report_only_the_omitted_data(self, flagged):
        summary = summarize_flags(flagged, flags=[DroppingComment.UNDEFINED_STEREOCHEMISTRY.value])

        assert list(summary.loc[summary["kind"] == "flag", "flag"]) == ["Undefined Stereochemistry"]
        assert count_of(summary, ANY_SELECTED_LABEL) == 2
        assert count_of(summary, RETAINED_LABEL) == 2

    def test_selected_flag_absent_from_data_is_reported_as_zero(self, flagged):
        summary = summarize_flags(flagged, flags=[DroppingComment.MIXTURE_IN_SMILES.value])

        assert count_of(summary, "Mixture in SMILES") == 0
        assert count_of(summary, ANY_SELECTED_LABEL) == 0
        assert count_of(summary, RETAINED_LABEL) == len(flagged)

    def test_thresholds_of_the_same_flag_are_pooled(self):
        data = pd.DataFrame({"data_dropping_comment": ["Assay size < 5", "Assay size < 20"]})

        summary = summarize_flags(data)

        assert count_of(summary, "Assay size <") == 2

    def test_selected_flags_are_normalized_before_matching(self):
        data = pd.DataFrame({"data_dropping_comment": ["Assay size < 5"]})

        summary = summarize_flags(data, flags=["Assay size < 5"])

        assert count_of(summary, "Assay size <") == 1

    def test_overlapping_flag_names_are_not_double_counted(self):
        """"Insufficient assay overlap" is a prefix of the metadata-matching variant."""
        data = pd.DataFrame(
            {
                "data_dropping_comment": [
                    "Insufficient assay overlap (min_overlap=5)",
                    "Insufficient assay overlap with metadata matching (min_overlap=5)",
                ]
            }
        )

        summary = summarize_flags(data)

        assert count_of(summary, "Insufficient assay overlap") == 1
        assert count_of(summary, "Insufficient assay overlap with metadata matching") == 1
        assert count_of(summary, ANY_FLAG_LABEL) == 2

    def test_flags_are_ordered_by_descending_count(self):
        data = pd.DataFrame(
            {"data_dropping_comment": ["Potential Duplicate"] * 3 + ["Undefined Stereochemistry"]}
        )

        summary = summarize_flags(data)

        assert list(summary.loc[summary["kind"] == "flag", "flag"]) == [
            "Potential Duplicate",
            "Undefined Stereochemistry",
        ]

    def test_empty_dataset_reports_zeros_without_dividing_by_zero(self):
        summary = summarize_flags(pd.DataFrame({"data_dropping_comment": []}))

        assert count_of(summary, TOTAL_LABEL) == 0
        assert pct_of(summary, ANY_FLAG_LABEL) == 0.0

    def test_missing_column_is_an_error(self):
        with pytest.raises(ValueError, match="not found in DataFrame"):
            summarize_flags(pd.DataFrame({"other": [1]}))

    def test_alternative_comment_column(self):
        """Pairwise comparisons carry their flags in `dropping_comment`."""
        data = pd.DataFrame({"dropping_comment": ["Undefined Stereochemistry|"]})

        summary = summarize_flags(data, comment_column="dropping_comment")

        assert count_of(summary, "Undefined Stereochemistry") == 1


class TestAggregatedCounting:
    @pytest.fixture
    def aggregated(self):
        """Two rows pooling five measurements, two of which are flagged."""
        return pd.DataFrame(
            {
                "data_dropping_comment": [
                    "Undefined Stereochemistry||Potential Duplicate",
                    "|",
                ]
            }
        )

    def test_row_counting_marks_a_row_carrying_any_flagged_measurement(self, aggregated):
        summary = summarize_flags(aggregated)

        assert count_of(summary, TOTAL_LABEL) == 2
        assert count_of(summary, ANY_FLAG_LABEL) == 1
        assert count_of(summary, "Undefined Stereochemistry") == 1

    def test_measurement_counting_uses_pooled_measurements_as_denominator(self, aggregated):
        summary = summarize_flags(aggregated, per_measurement=True)

        assert count_of(summary, TOTAL_LABEL) == 5
        assert count_of(summary, ANY_FLAG_LABEL) == 2
        assert count_of(summary, UNFLAGGED_LABEL) == 3
        assert pct_of(summary, ANY_FLAG_LABEL) == 40.0

    def test_row_and_measurement_counts_differ_on_aggregated_data(self, aggregated):
        by_row = summarize_flags(aggregated)
        by_measurement = summarize_flags(aggregated, per_measurement=True)

        assert count_of(by_row, TOTAL_LABEL) < count_of(by_measurement, TOTAL_LABEL)


class TestSummaryMatchesFiltering:
    """The reported omission must equal what dropping those flags actually removes."""

    def test_union_equals_rows_removed_by_filter_dropping_flags(self):
        from Capricho.core.pandas_helper import filter_dropping_flags

        data = pd.DataFrame(
            {
                "data_dropping_comment": [
                    "Undefined Stereochemistry",
                    "Undefined Stereochemistry & Potential Duplicate",
                    "Potential Duplicate",
                    "",
                    "Mixture in SMILES",
                ]
            }
        )
        flags = [
            DroppingComment.UNDEFINED_STEREOCHEMISTRY.value,
            DroppingComment.POTENTIAL_DUPLICATE.value,
        ]

        summary = summarize_flags(data, flags=flags)
        filtered = filter_dropping_flags(data, flags)

        assert count_of(summary, ANY_SELECTED_LABEL) == len(data) - len(filtered)
        assert count_of(summary, RETAINED_LABEL) == len(filtered)


class TestFormatFlagSummary:
    def test_reports_absolute_number_and_percentage(self):
        data = pd.DataFrame({"data_dropping_comment": ["Potential Duplicate", ""]})

        text = format_flag_summary(summarize_flags(data))

        assert "Potential Duplicate:" in text
        assert "1" in text and "50.0%" in text

    def test_totals_are_separated_from_the_per_flag_rows(self):
        data = pd.DataFrame({"data_dropping_comment": ["Potential Duplicate"]})

        lines = format_flag_summary(summarize_flags(data)).splitlines()

        assert lines[0].strip() == "QUALITY FLAGS"
        assert any(set(line.strip()) == {"-"} for line in lines)
        assert ANY_FLAG_LABEL in lines[-3]
        assert TOTAL_LABEL in lines[-1]

    def test_unflagged_dataset_says_so(self):
        data = pd.DataFrame({"data_dropping_comment": ["", ""]})

        text = format_flag_summary(summarize_flags(data))

        assert "(none)" in text


class TestCrossAssayCoverage:
    @pytest.fixture
    def aggregated(self):
        """Three compounds; only the first two were measured in more than one assay."""
        return pd.DataFrame(
            {"assay_chembl_id": ["CHEMBL1|CHEMBL2|CHEMBL3", "CHEMBL1|CHEMBL2", "CHEMBL1"]}
        )

    def test_counts_datapoints_reachable_by_cross_assay_comparison(self, aggregated):
        coverage = summarize_cross_assay_coverage(aggregated)

        assert coverage["aggregated_datapoints"] == 3
        assert coverage["comparable_datapoints"] == 2
        assert coverage["pct_comparable_datapoints"] == pytest.approx(66.67, abs=0.01)

    def test_counts_the_measurements_those_datapoints_pool(self, aggregated):
        coverage = summarize_cross_assay_coverage(aggregated)

        assert coverage["measurements_in_comparable"] == 5

    def test_share_of_retrieved_measurements_when_known(self, aggregated):
        coverage = summarize_cross_assay_coverage(aggregated, n_retrieved=10)

        assert coverage["pct_measurements_in_comparable"] == 50.0

    def test_share_of_retrieved_is_omitted_when_unknown(self, aggregated):
        assert "pct_measurements_in_comparable" not in summarize_cross_assay_coverage(aggregated)

    def test_single_assay_dataset_has_no_coverage(self):
        coverage = summarize_cross_assay_coverage(pd.DataFrame({"assay_chembl_id": ["CHEMBL1"]}))

        assert coverage["comparable_datapoints"] == 0
        assert coverage["pct_comparable_datapoints"] == 0.0

    def test_missing_column_is_an_error(self):
        with pytest.raises(ValueError, match="not found in DataFrame"):
            summarize_cross_assay_coverage(pd.DataFrame({"other": [1]}))

    def test_formats_counts_and_shares(self, aggregated):
        text = format_cross_assay_coverage(summarize_cross_assay_coverage(aggregated, n_retrieved=10))

        assert "Measured in >1 assay:" in text
        assert "66.7%" in text
        assert "50.0% of retrieved" in text
