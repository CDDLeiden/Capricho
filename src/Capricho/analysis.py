"""Analysis utilities for CAPRICHO bioactivity data comparability studies"""

from enum import Enum
from itertools import combinations
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from scipy import stats
from sklearn.metrics import r2_score


class ProcessingComment(str, Enum):
    """Processing comments added during data curation (non-dropping flags)."""

    CALCULATED_PCHEMBL = "Calculated pChEMBL"
    SALT_SOLVENT_REMOVED = "Salt/solvent removed"
    PCHEMBL_DUPLICATION_ACROSS_DOCUMENTS = "pChEMBL Duplication Across Documents"


class DroppingComment(str, Enum):
    """Dropping comments indicating data quality issues (dropping flags).

    Note: Assay size values shown are examples. Actual thresholds depend on user
    parameters and should be matched using pattern-based functions.
    """

    DATA_VALIDITY_COMMENT = "Data Validity Comment Present"
    POTENTIAL_DUPLICATE = "Potential Duplicate"
    UNDEFINED_STEREOCHEMISTRY = "Undefined Stereochemistry"
    MUTATION_KEYWORD = "Mutation keyword in assay description"
    ASSAY_SIZE_TOO_LARGE = "Assay size >"  # Example: "Assay size > 100"
    ASSAY_SIZE_TOO_SMALL = "Assay size <"  # Example: "Assay size < 20"
    UNIT_ANNOTATION_ERROR = "Unit Annotation Error"


def normalize_comment_pattern(comment: str) -> str:
    """Normalize comment by removing dynamic values (e.g., threshold numbers).

    Args:
        comment: Raw comment string from data.

    Returns:
        Normalized pattern for matching. For example:
        - "Assay size < 20" -> "Assay size <"
        - "Assay size > 100" -> "Assay size >"
        - "Unit Annotation Error" -> "Unit Annotation Error"
    """
    # Handle assay size patterns by removing the threshold number
    if comment.startswith("Assay size <"):
        return "Assay size <"
    elif comment.startswith("Assay size >"):
        return "Assay size >"
    return comment


def extract_assay_threshold(comment: str) -> str:
    """Extract the threshold value from an assay size comment.

    Args:
        comment: Comment like "Assay size < 20" or "Assay size > 100".

    Returns:
        The threshold number as string, or empty string if not found.
    """
    if comment.startswith("Assay size <"):
        return comment.split("<")[1].strip()
    elif comment.startswith("Assay size >"):
        return comment.split(">")[1].strip()
    return ""


def get_all_comments() -> list[str]:
    """Get all comment types in order for plotting.

    Returns:
        List of all comment patterns (dropping and processing combined).
        Assay size comments are patterns without specific thresholds.
    """
    return [
        DroppingComment.DATA_VALIDITY_COMMENT.value,
        DroppingComment.POTENTIAL_DUPLICATE.value,
        DroppingComment.UNDEFINED_STEREOCHEMISTRY.value,
        DroppingComment.MUTATION_KEYWORD.value,
        DroppingComment.ASSAY_SIZE_TOO_LARGE.value,
        DroppingComment.ASSAY_SIZE_TOO_SMALL.value,
        DroppingComment.UNIT_ANNOTATION_ERROR.value,
        ProcessingComment.SALT_SOLVENT_REMOVED.value,
        ProcessingComment.CALCULATED_PCHEMBL.value,
        ProcessingComment.PCHEMBL_DUPLICATION_ACROSS_DOCUMENTS.value,
    ]


def deaggregate_data(data: pd.DataFrame, sep_str: str = "|") -> pd.DataFrame:
    """De-aggregate data by splitting pipe-delimited columns back into individual rows.

    Takes aggregated data where multiple values are stored in pipe-delimited strings
    and explodes them back into separate rows, effectively reversing the aggregation.
    Only rows containing the separator in at least one column are exploded.

    Args:
        data: DataFrame with potentially aggregated (pipe-delimited) columns.
        sep_str: Separator string used to delimit multiple values in columns.

    Returns:
        DataFrame with de-aggregated data where each measurement is a separate row.
    """
    cols_with_pipe = data.apply(lambda col: col.astype(str).str.contains(sep_str, regex=False)).any()
    cols_with_pipe = np.compress(cols_with_pipe.values, cols_with_pipe.index).tolist()

    if not cols_with_pipe:
        return data.copy()

    first_pipe_col = cols_with_pipe[0]
    mask = data[first_pipe_col].astype(str).str.contains(sep_str, regex=False)
    aggregated_rows = data[mask]

    if len(aggregated_rows) == 0:
        return data.copy()

    aggregated_rows = aggregated_rows.copy()
    aggregated_rows.loc[:, cols_with_pipe] = aggregated_rows.loc[:, cols_with_pipe].apply(
        lambda col: col.fillna("").astype(str).str.split(sep_str)
    )

    exploded = aggregated_rows.explode(column=cols_with_pipe)
    deaggregated = pd.concat([data.drop(index=aggregated_rows.index), exploded])

    return deaggregated.reset_index(drop=True)


def explode_assay_comparability(subset: pd.DataFrame, sep_str: str = "|") -> pd.DataFrame:
    """Explode dataset to create pairwise comparisons between assays for the same compound.

    Takes a subset of data where compounds have measurements across multiple assays (indicated
    by separator in columns) and creates all pairwise combinations for comparability analysis.

    Args:
        subset: DataFrame with multi-valued columns separated by sep_str.
        sep_str: Separator string used to delimit multiple values in columns.

    Returns:
        DataFrame with exploded pairwise comparisons, with _x and _y suffixes for each pair.
    """
    singleval_cols = [
        "connectivity",
        "target_chembl_id",
        "repeat",
    ]
    multival_cols = [
        "activity_id",
        "assay_chembl_id",
        "pchembl_value",
        "data_processing_comment",
        "data_dropping_comment",
        "standard_type",
    ]

    exploded_subset = subset[
        [
            *singleval_cols,
            *multival_cols,
        ]
    ].apply(lambda x: x.str.split(sep_str) if x.name not in singleval_cols else x)

    for col in multival_cols:
        exploded_subset[col] = exploded_subset[col].apply(
            lambda x: [sep_str.join(y) for y in combinations(x, 2)]
        )

    exploded_subset = exploded_subset.explode(multival_cols)

    # Handle empty DataFrame case
    if len(exploded_subset) == 0:
        # Create empty result with correct columns
        result_cols = [*singleval_cols]
        for col in multival_cols:
            result_cols.extend([f"{col}_x", f"{col}_y"])
        result_cols.extend(["processing_comment", "dropping_comment"])
        return pd.DataFrame(columns=result_cols)

    suffixes = ["_x", "_y"]
    for col in multival_cols:
        values = exploded_subset[col].apply(lambda x: x.split(sep_str)).values
        for idx, s in enumerate(suffixes):
            exploded_subset[f"{col}{s}"] = [v[idx] for v in values]
        exploded_subset.drop(columns=[col], inplace=True)

    # Combine processing and dropping comments from both assays
    exploded_subset = exploded_subset.assign(
        processing_comment=lambda x: x["data_processing_comment_x"].fillna("")
        + sep_str
        + x["data_processing_comment_y"].fillna(""),
        dropping_comment=lambda x: x["data_dropping_comment_x"].fillna("")
        + sep_str
        + x["data_dropping_comment_y"].fillna(""),
    ).query("assay_chembl_id_x != assay_chembl_id_y")

    # Clean up comment strings by removing duplicates and extra separators
    exploded_subset["dropping_comment"] = (
        exploded_subset.dropping_comment.str.rstrip(sep_str)
        .str.lstrip(sep_str)
        .apply(lambda x: sep_str.join(np.unique(x.split(sep_str))))
    )
    exploded_subset["processing_comment"] = (
        exploded_subset.processing_comment.str.rstrip(sep_str)
        .str.lstrip(sep_str)
        .apply(lambda x: sep_str.join(np.unique(x.split(sep_str))))
    )
    return exploded_subset


def plot_subset(
    subset: pd.DataFrame,
    title: str = "",
    color: str = "slategray",
    figsize: Tuple[float, float] = (5, 5),
) -> Tuple[plt.Figure, plt.Axes]:
    """Create scatter plot comparing pChEMBL values across assays with correlation metrics.

    Args:
        subset: DataFrame with pchembl_value_x and pchembl_value_y columns.
        title: Plot title.
        color: Color for scatter points.
        figsize: Figure size as (width, height) tuple.

    Returns:
        Tuple of (figure, axes) objects.
    """
    fig, ax = plt.subplots(figsize=figsize)
    subset = subset.copy()
    subset["pchembl_value_x"] = subset["pchembl_value_x"].astype(float)
    subset["pchembl_value_y"] = subset["pchembl_value_y"].astype(float)

    ax.scatter(
        subset["pchembl_value_x"],
        subset["pchembl_value_y"],
        alpha=0.3,
        edgecolors="none",
        color=color,
    )
    ax.set_title(title)

    # Add reference lines
    ax.plot((3, 12), (3, 12), "k-", label="Identity $(y=x)$")
    ax.plot((3, 12), (2, 11), "k--")
    ax.plot((3, 12), (4, 13), "k--", label="$y=x±1$")
    ax.plot((3, 12), (2.7, 11.7), "k-.")
    ax.plot((3, 12), (3.3, 12.3), "k-.", label="$y=x±0.3$")

    ax.set_xlim(3, 12)
    ax.set_ylim(3, 12)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Calculate correlation metrics
    xp = subset["pchembl_value_x"]
    yp = subset["pchembl_value_y"]

    r, _ = stats.spearmanr(xp, yp)
    tau, _ = stats.kendalltau(xp, yp)
    r2 = r2_score(xp, yp)

    ax.text(
        1.0,
        0.175,
        rf"$R^2: {r2:.2f}$" + "\n" + rf"Spearman $\rho: {r:.2f}$" + "\n" + rf"Kendall $\tau: {tau:.2f}$",
        transform=ax.transAxes,
        verticalalignment="top",
        horizontalalignment="right",
    )

    ax.set_ylabel("Assay 2 pChEMBL value")
    ax.set_xlabel("Assay 1 pChEMBL value")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[1:],
        labels[1:],
        title="",
        bbox_to_anchor=(1.05, 0, 1, 0.2),
        loc="lower left",
        mode="expand",
        frameon=False,
    )
    fig.tight_layout()

    return fig, ax


def build_query_string(comment: str) -> str:
    """Build query string for filtering data by specific comment flags.

    This function replicates the original notebook logic for querying exploded datasets.
    It handles both exact matches and pattern-based matches (e.g., for assay size).

    Args:
        comment: The comment pattern to search for (may be normalized pattern like "Assay size <").

    Returns:
        Query string suitable for pd.DataFrame.query().
    """
    # Special handling for document duplication (needs to be in both assays)
    if comment == ProcessingComment.PCHEMBL_DUPLICATION_ACROSS_DOCUMENTS.value:
        return (
            "(data_processing_comment_x.str.contains('pChEMBL Duplication Across Documents', regex=False) & "
            "data_processing_comment_y.str.contains('pChEMBL Duplication Across Documents', regex=False)) | "
            "(data_processing_comment_y.str.contains('pChEMBL Duplication Across Documents', regex=False) & "
            "data_processing_comment_x.str.contains('pChEMBL Duplication Across Documents', regex=False))"
        )

    # Special handling for calculated pChEMBL (uses combined processing_comment column)
    if comment == ProcessingComment.CALCULATED_PCHEMBL.value:
        return "processing_comment.str.contains('Calculated pChEMBL', regex=False) & (dropping_comment == '')"

    # For other processing comments (uses combined processing_comment column)
    if comment in [ProcessingComment.SALT_SOLVENT_REMOVED.value]:
        return f"processing_comment.str.contains('{comment}', regex=False) & dropping_comment == ''"

    # Special handling for potential duplicates (needs to be in both assays, exclude data validity issues)
    if comment == DroppingComment.POTENTIAL_DUPLICATE.value:
        return (
            "((data_dropping_comment_x.str.contains('Potential Duplicate', regex=False) & "
            "data_dropping_comment_y.str.contains('Potential Duplicate', regex=False)) | "
            "(data_dropping_comment_y.str.contains('Potential Duplicate', regex=False) & "
            "data_dropping_comment_x.str.contains('Potential Duplicate', regex=False))) & "
            "~dropping_comment.str.contains('Data Validity Comment Present', regex=False)"
        )

    # Special handling for unit annotation errors (needs to be in both assays, exclude data validity issues)
    if comment == DroppingComment.UNIT_ANNOTATION_ERROR.value:
        return (
            "((data_dropping_comment_x.str.contains('Unit Annotation Error', regex=False) & "
            "data_dropping_comment_y.str.contains('Unit Annotation Error', regex=False)) | "
            "(data_dropping_comment_y.str.contains('Unit Annotation Error', regex=False) & "
            "data_dropping_comment_x.str.contains('Unit Annotation Error', regex=False))) & "
            "~dropping_comment.str.contains('Data Validity Comment Present', regex=False)"
        )

    # For Data Validity Comment, allow any processing comment
    if comment == DroppingComment.DATA_VALIDITY_COMMENT.value:
        return f"dropping_comment.str.contains('{comment}', regex=False)"

    # For other dropping comments, exclude processing comments and Data Validity Comment
    return (
        f"dropping_comment.str.contains('{comment}', regex=False) & "
        "processing_comment == '' & "
        "~dropping_comment.str.contains('Data Validity Comment Present', regex=False)"
    )


def plot_multi_panel_comparability(
    exploded_subset: pd.DataFrame,
    comments: List[str],
    title: str = "Comparability Across Flagged Data",
    figsize: Tuple[float, float] = (20, 8),
    ncols: int = 5,
) -> Tuple[plt.Figure, np.ndarray]:
    """Create multi-panel plot showing comparability for different data quality flags.

    Args:
        exploded_subset: DataFrame from explode_assay_comparability().
        comments: List of comment strings to plot.
        title: Overall figure title.
        figsize: Figure size as (width, height) tuple.
        ncols: Number of columns in subplot grid.

    Returns:
        Tuple of (figure, axes array).
    """
    nrows = int(np.ceil(len(comments) / ncols))
    colors = [tuple([*col] + [1]) for col in colormaps["tab10"].colors]

    fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
    axs_flat = axs.flatten()

    for idx, color, obs, ax in zip(range(1, len(comments) + 1), colors, comments, axs_flat):
        query_str = build_query_string(obs)
        subset = exploded_subset.query(query_str)

        # Extract actual title with dynamic threshold if it's an assay size comment
        title_str = obs
        if obs.startswith("Assay size"):
            # Find first occurrence with actual threshold from data
            if len(subset) > 0:
                for col in ["data_dropping_comment_x", "data_dropping_comment_y"]:
                    if col in subset.columns:
                        sample_comments = subset[col].dropna()
                        if len(sample_comments) > 0:
                            for comment in sample_comments.iloc[0].split("|"):
                                if normalize_comment_pattern(comment) == obs:
                                    title_str = comment
                                    break
                            if title_str != obs:
                                break

        if len(subset) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(f"{idx}. {title_str}")
            continue

        subset = subset.copy()
        subset["pchembl_value_x"] = subset["pchembl_value_x"].astype(float)
        subset["pchembl_value_y"] = subset["pchembl_value_y"].astype(float)

        ax.scatter(
            subset["pchembl_value_x"],
            subset["pchembl_value_y"],
            alpha=0.3,
            edgecolors="none",
            label=title_str,
            color=color,
        )
        ax.set_title(f"{idx}. {title_str}")

        # Add reference lines
        ax.plot((3, 12), (3, 12), "k-", label="Identity $(y=x)$")
        ax.plot((3, 12), (2, 11), "k--")
        ax.plot((3, 12), (4, 13), "k--", label="$y=x±1$")
        ax.plot((3, 12), (2.7, 11.7), "k-.")
        ax.plot((3, 12), (3.3, 12.3), "k-.", label="$y=x±0.3$")

        ax.set_xlim(3, 12)
        ax.set_ylim(3, 12)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        xp = subset["pchembl_value_x"]
        yp = subset["pchembl_value_y"]

        try:
            r, _ = stats.spearmanr(xp, yp)
            tau, _ = stats.kendalltau(xp, yp)
            r2 = r2_score(xp, yp)
        except ValueError:
            continue

        ax.text(
            1.0,
            0.175,
            rf"$R^2: {r2:.2f}$" + "\n" + rf"Spearman $\rho: {r:.2f}$" + "\n" + rf"Kendall $\tau: {tau:.2f}$",
            transform=ax.transAxes,
            verticalalignment="top",
            horizontalalignment="right",
        )

        if idx in [1, ncols + 1]:
            ax.set_ylabel("Assay 2 pChEMBL value")
        if idx > (nrows - 1) * ncols:
            ax.set_xlabel("Assay 1 pChEMBL value")

        if idx == len(comments):
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(
                handles[1:],
                labels[1:],
                title="",
                bbox_to_anchor=(1.05, 0, 1, 0.2),
                loc="lower left",
                mode="expand",
                frameon=False,
            )

    # Hide unused subplots
    for ax in axs_flat[len(comments) :]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.suptitle(title, fontsize=16, y=1.02)

    return fig, axs
