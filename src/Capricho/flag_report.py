"""Account for how much data carries each quality flag and how much is omitted when flags are dropped.

CAPRICHO flags data rather than silently discarding it, so the amount of data a user
actually loses depends on which flags they choose to drop. These helpers separate the
two questions:

- *flagged*: how many entries carry each flag, reported per flag for transparency;
- *omitted*: how many entries carry at least one of the flags being dropped.

The second is not the sum of the first. Flags co-occur — a measurement can be flagged
for a small assay and undefined stereochemistry at once — so summing per-flag counts
overcounts the loss. The ``Any flag`` total is the union and is the number to quote as
"data omitted".

Percentages are always relative to the entries in the DataFrame that is passed in, so
the caller decides the denominator: retrieved measurements, aggregated datapoints, or
pairwise cross-assay comparisons.
"""

from typing import Optional, Union

import pandas as pd

from .core.default_fields import DATA_DROPPING_COMMENT

#: Separator joining several flags written onto the same measurement (see ``add_comment``).
FLAG_SEPARATOR = " & "

ANY_FLAG_LABEL = "Any flag"
ANY_SELECTED_LABEL = "Any of the selected flags"
UNFLAGGED_LABEL = "Unflagged"
RETAINED_LABEL = "Retained"
TOTAL_LABEL = "All entries"


def measurement_comments(cell, sep_str: str = "|") -> list[str]:
    """Split one comment cell into the comment of each measurement it pools.

    Aggregated rows hold one comment per pooled measurement, separated by ``sep_str``.
    Non-aggregated rows yield a single comment.

    Args:
        cell: Value of a comment column for one row. Missing values yield one empty comment.
        sep_str: Separator delimiting pooled measurements.

    Returns:
        List of per-measurement comment strings, always at least one element long.
    """
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return [""]
    text = str(cell)
    if text == "nan":
        return [""]
    return text.split(sep_str)


def flags_in_comment(comment: str) -> set[str]:
    """Extract the normalized flag names present in a single measurement's comment.

    Flags carrying a run-specific value are normalized to their pattern, so that
    "Assay size < 5" and "Assay size < 20" are counted as the same reason.

    Args:
        comment: Comment of a single measurement, possibly joining several flags.

    Returns:
        Set of normalized flag names. Empty for an empty or missing comment.
    """
    from .analysis import normalize_comment_pattern

    if not comment or comment == "nan":
        return set()
    return {
        normalize_comment_pattern(flag.strip()) for flag in comment.split(FLAG_SEPARATOR) if flag.strip()
    }


def _entry_flags(
    data: pd.DataFrame,
    comment_column: str,
    sep_str: str,
    per_measurement: bool,
) -> list[set[str]]:
    """Reduce a comment column to one set of normalized flags per counted entry."""
    entries: list[set[str]] = []
    for cell in data[comment_column]:
        comments = measurement_comments(cell, sep_str)
        if per_measurement:
            entries.extend(flags_in_comment(comment) for comment in comments)
        else:
            row_flags: set[str] = set()
            for comment in comments:
                row_flags |= flags_in_comment(comment)
            entries.append(row_flags)
    return entries


def summarize_flags(
    data: pd.DataFrame,
    flags: Optional[list[str]] = None,
    comment_column: str = DATA_DROPPING_COMMENT,
    sep_str: str = "|",
    per_measurement: bool = False,
) -> pd.DataFrame:
    """Count how many entries carry each quality flag, in absolute numbers and as a share.

    Args:
        data: DataFrame holding a comment column. May be non-aggregated (one row per
            measurement), aggregated (comments pooled with ``sep_str``), or exploded
            into pairwise comparisons.
        flags: Flags to account for, normally the ones being dropped. Their union is
            reported as the data omitted. If None, every flag found in the data is
            reported and the union covers all of them.
        comment_column: Column holding the flags. Defaults to ``data_dropping_comment``;
            pass ``"dropping_comment"`` for data exploded by
            :func:`~Capricho.analysis.explode_assay_comparability`.
        sep_str: Separator delimiting pooled measurements within an aggregated cell.
        per_measurement: Count each pooled measurement separately instead of counting
            rows. Meaningful for aggregated data, where one row can pool many
            measurements. Leave False for pairwise comparisons, whose unit is the pair.

    Returns:
        DataFrame with columns ``flag``, ``kind``, ``n`` and ``pct``, sorted by
        descending count. Rows with ``kind == "flag"`` are the per-flag counts; the
        three ``kind == "total"`` rows are the union of the reported flags, the entries
        carrying none of them, and the denominator. The union is labelled
        ``"Any of the selected flags"`` and its complement ``"Retained"`` when ``flags``
        is given, describing a filter; without a selection they are ``"Any flag"`` and
        ``"Unflagged"``, describing the data as it stands.

    Raises:
        ValueError: If ``comment_column`` is not a column of ``data``.
    """
    from .analysis import normalize_comment_pattern

    if comment_column not in data.columns:
        raise ValueError(f"Column '{comment_column}' not found in DataFrame. Columns: {list(data.columns)}")

    entries = _entry_flags(data, comment_column, sep_str, per_measurement)
    n_total = len(entries)

    if flags is None:
        names = sorted({flag for entry in entries for flag in entry})
        union_label, complement_label = ANY_FLAG_LABEL, UNFLAGGED_LABEL
    else:
        names = list(dict.fromkeys(normalize_comment_pattern(str(flag)) for flag in flags))
        union_label, complement_label = ANY_SELECTED_LABEL, RETAINED_LABEL
    selected = set(names)

    counts = {name: sum(1 for entry in entries if name in entry) for name in names}
    n_omitted = sum(1 for entry in entries if entry & selected)

    def _pct(count: int) -> float:
        return count / n_total * 100 if n_total else 0.0

    rows = [
        {"flag": name, "kind": "flag", "n": counts[name], "pct": _pct(counts[name])}
        for name in sorted(names, key=lambda name: (-counts[name], name))
    ]
    rows += [
        {"flag": union_label, "kind": "total", "n": n_omitted, "pct": _pct(n_omitted)},
        {
            "flag": complement_label,
            "kind": "total",
            "n": n_total - n_omitted,
            "pct": _pct(n_total - n_omitted),
        },
        {"flag": TOTAL_LABEL, "kind": "total", "n": n_total, "pct": 100.0 if n_total else 0.0},
    ]
    return pd.DataFrame(rows, columns=["flag", "kind", "n", "pct"])


def format_flag_summary(
    summary: pd.DataFrame,
    title: str = "QUALITY FLAGS",
    indent: str = "  ",
    empty_message: str = "(none)",
) -> str:
    """Render a flag summary as an aligned text block for logging or printing.

    Args:
        summary: Output of :func:`summarize_flags`.
        title: Heading placed above the table.
        indent: Prefix for the heading; table rows are indented one level deeper.
        empty_message: Line shown in place of the table when no flag was found.

    Returns:
        Multi-line string, without a trailing newline.
    """
    per_flag = summary[summary["kind"] == "flag"]
    totals = summary[summary["kind"] == "total"]

    lines = [f"{indent}{title}"]
    row_indent = indent * 2
    width = max(max((len(str(flag)) for flag in summary["flag"]), default=0) + 1, 30)

    def _row(flag: str, n: int, pct: float) -> str:
        return f"{row_indent}{flag + ':':<{width}s} {n:>8,}  ({pct:5.1f}%)"

    if per_flag.empty:
        lines.append(f"{row_indent}{empty_message}")
    for row in per_flag.itertuples():
        lines.append(_row(row.flag, row.n, row.pct))
    lines.append(f"{row_indent}{'-' * (width + 20)}")
    for row in totals.itertuples():
        lines.append(_row(row.flag, row.n, row.pct))
    return "\n".join(lines)


def summarize_cross_assay_coverage(
    data: pd.DataFrame,
    assay_column: str = "assay_chembl_id",
    sep_str: str = "|",
    n_retrieved: Optional[int] = None,
) -> dict[str, Union[int, float]]:
    """Measure how much of an aggregated dataset cross-assay comparability can reach.

    Only compounds measured in more than one assay contribute a comparison, so this is
    the share of the data the comparability analysis actually observes.

    Args:
        data: Aggregated DataFrame, one row per compound-target readout.
        assay_column: Column holding the pooled assay identifiers.
        sep_str: Separator delimiting pooled measurements.
        n_retrieved: Number of measurements retrieved before aggregation, if known.
            Adds the share of retrieved measurements that the comparable rows pool.

    Returns:
        Dict with the number of aggregated datapoints, how many are measured in more
        than one assay, that count as a percentage, the number of measurements those
        rows pool, and — when ``n_retrieved`` is given — that count as a percentage of
        the retrieved measurements.

    Raises:
        ValueError: If ``assay_column`` is not a column of ``data``.
    """
    if assay_column not in data.columns:
        raise ValueError(f"Column '{assay_column}' not found in DataFrame. Columns: {list(data.columns)}")

    n_assays = data[assay_column].map(lambda cell: len(measurement_comments(cell, sep_str)))
    multi_assay = n_assays > 1
    n_datapoints = len(data)
    n_comparable = int(multi_assay.sum())

    coverage = {
        "aggregated_datapoints": n_datapoints,
        "comparable_datapoints": n_comparable,
        "pct_comparable_datapoints": (n_comparable / n_datapoints * 100 if n_datapoints else 0.0),
        "measurements_in_comparable": int(n_assays[multi_assay].sum()),
    }
    if n_retrieved:
        coverage["pct_measurements_in_comparable"] = (
            coverage["measurements_in_comparable"] / n_retrieved * 100
        )
    return coverage


def format_cross_assay_coverage(
    coverage: dict[str, Union[int, float]],
    title: str = "CROSS-ASSAY COVERAGE",
    indent: str = "  ",
) -> str:
    """Render :func:`summarize_cross_assay_coverage` output as an aligned text block.

    Args:
        coverage: Output of :func:`summarize_cross_assay_coverage`.
        title: Heading placed above the block.
        indent: Prefix for the heading; entries are indented one level deeper.

    Returns:
        Multi-line string, without a trailing newline.
    """
    row_indent = indent * 2
    width = 32
    lines = [
        f"{indent}{title}",
        f"{row_indent}{'Aggregated datapoints:':<{width}s} {coverage['aggregated_datapoints']:>8,}",
        f"{row_indent}{'Measured in >1 assay:':<{width}s} {coverage['comparable_datapoints']:>8,}"
        f"  ({coverage['pct_comparable_datapoints']:5.1f}%)",
    ]
    line = f"{row_indent}{'Measurements they pool:':<{width}s} {coverage['measurements_in_comparable']:>8,}"
    if "pct_measurements_in_comparable" in coverage:
        line += f"  ({coverage['pct_measurements_in_comparable']:5.1f}% of retrieved)"
    lines.append(line)
    return "\n".join(lines)
