import marimo

__generated_with = "0.20.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import io

    import altair as alt
    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    from scipy import stats
    from sklearn.metrics import r2_score

    from Capricho.analysis import (
        DroppingComment,
        ProcessingComment,
        build_query_string,
        explode_assay_comparability,
        get_all_comments,
        plot_multi_panel_comparability,
    )
    from Capricho.cli.prepare import clean_data
    from Capricho.core.binarization import (
        VALID_CONFLICT_STRATEGIES,
        binarize_aggregated_data,
    )
    from Capricho.core.default_fields import (
        DATA_DROPPING_COMMENT,
        DATA_PROCESSING_COMMENT,
        multiple_value_cols,
    )

    return (
        AllChem,
        Chem,
        DATA_DROPPING_COMMENT,
        DATA_PROCESSING_COMMENT,
        DroppingComment,
        VALID_CONFLICT_STRATEGIES,
        alt,
        binarize_aggregated_data,
        build_query_string,
        clean_data,
        explode_assay_comparability,
        get_all_comments,
        io,
        mo,
        multiple_value_cols,
        np,
        pd,
        plot_multi_panel_comparability,
        plt,
        r2_score,
        rdMolDraw2D,
        stats,
    )


@app.cell
def _(
    AllChem,
    Chem,
    DATA_DROPPING_COMMENT,
    DATA_PROCESSING_COMMENT,
    alt,
    get_all_comments,
    mo,
    multiple_value_cols,
    pd,
    r2_score,
    rdMolDraw2D,
    stats,
):
    def smiles_to_svg_html(smiles):
        if pd.isna(smiles) or smiles == "":
            return ""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""

        AllChem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DSVG(150, 150)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        svg = drawer.GetDrawingText()
        return mo.Html(svg)

    def detect_data_properties(df):
        """Analyze a loaded DataFrame to detect its structure and contents."""
        is_aggregated = False
        for col in multiple_value_cols:
            if col in df.columns and df[col].astype(str).str.contains("|", regex=False).any():
                is_aggregated = True
                break

        value_column = "pchembl_value"
        if "pchembl_value" not in df.columns:
            if "standard_value" in df.columns:
                value_column = "standard_value"

        mean_col = f"{value_column}_mean"
        if mean_col in df.columns:
            vals = pd.to_numeric(df[mean_col], errors="coerce").dropna()
            if len(vals) > 0:
                margin = max((vals.max() - vals.min()) * 0.1, 0.5)
                axis_limits = (float(vals.min() - margin), float(vals.max() + margin))
            else:
                axis_limits = (3, 12) if value_column == "pchembl_value" else None
        elif value_column == "pchembl_value":
            axis_limits = (3, 12)
        else:
            axis_limits = None

        present_flags = []
        all_flags = get_all_comments()
        for flag in all_flags:
            for col in [DATA_DROPPING_COMMENT, DATA_PROCESSING_COMMENT]:
                if col in df.columns and df[col].fillna("").str.contains(flag, regex=False).any():
                    present_flags.append(flag)
                    break

        return {
            "is_aggregated": is_aggregated,
            "value_column": value_column,
            "axis_limits": axis_limits,
            "present_flags": present_flags,
            "has_dropping_flags": DATA_DROPPING_COMMENT in df.columns,
            "has_processing_flags": DATA_PROCESSING_COMMENT in df.columns,
        }

    def load_dataframe(source, filename):
        """Load a DataFrame from bytes or a file path, detecting format from extension."""
        if filename.endswith(".parquet"):
            return pd.read_parquet(source)
        elif filename.endswith(".tsv"):
            return pd.read_csv(source, sep="\t")
        else:
            return pd.read_csv(source)

    def build_comparability_chart(chart_data, title, color, value_column, axis_limits):
        """Build an Altair scatter plot with identity and reference lines for comparability analysis."""
        x_col = f"{value_column}_x"
        y_col = f"{value_column}_y"
        lo, hi = axis_limits

        base = (
            alt.Chart(chart_data)
            .mark_circle(size=60, opacity=0.5, color=color)
            .encode(
                x=alt.X(f"{x_col}:Q", title="Assay 1 value", scale=alt.Scale(domain=[lo, hi])),
                y=alt.Y(f"{y_col}:Q", title="Assay 2 value", scale=alt.Scale(domain=[lo, hi])),
                tooltip=[
                    alt.Tooltip("connectivity:N", title="Compound"),
                    alt.Tooltip(f"{x_col}:Q", title="Assay 1", format=".2f"),
                    alt.Tooltip(f"{y_col}:Q", title="Assay 2", format=".2f"),
                    alt.Tooltip("assay_chembl_id_x:N", title="Assay 1 ID"),
                    alt.Tooltip("assay_chembl_id_y:N", title="Assay 2 ID"),
                ],
            )
        )

        identity = (
            alt.Chart(pd.DataFrame({"x": [lo, hi], "y": [lo, hi]}))
            .mark_line(color="black", strokeDash=[0])
            .encode(x="x:Q", y="y:Q")
        )

        plus_one = (
            alt.Chart(pd.DataFrame({"x": [lo, hi], "y": [lo + 1, hi + 1]}))
            .mark_line(color="black", strokeDash=[5, 5])
            .encode(x="x:Q", y="y:Q")
        )

        minus_one = (
            alt.Chart(pd.DataFrame({"x": [lo, hi], "y": [lo - 1, hi - 1]}))
            .mark_line(color="black", strokeDash=[5, 5])
            .encode(x="x:Q", y="y:Q")
        )

        plus_03 = (
            alt.Chart(pd.DataFrame({"x": [lo, hi], "y": [lo + 0.3, hi + 0.3]}))
            .mark_line(color="black", strokeDash=[2, 2])
            .encode(x="x:Q", y="y:Q")
        )

        minus_03 = (
            alt.Chart(pd.DataFrame({"x": [lo, hi], "y": [lo - 0.3, hi - 0.3]}))
            .mark_line(color="black", strokeDash=[2, 2])
            .encode(x="x:Q", y="y:Q")
        )

        return (base + identity + plus_one + minus_one + plus_03 + minus_03).properties(
            width=600, height=600, title=title
        )

    def compute_correlation_metrics(x, y):
        """Compute R², Spearman ρ, and Kendall τ for two Series."""
        r, _ = stats.spearmanr(x, y)
        tau, _ = stats.kendalltau(x, y)
        r2 = r2_score(x, y)
        return {"r2": r2, "spearman": r, "kendall": tau}

    return (
        build_comparability_chart,
        compute_correlation_metrics,
        detect_data_properties,
        load_dataframe,
        smiles_to_svg_html,
    )


@app.cell
def _(mo):
    file_upload = mo.ui.file(
        filetypes=[".csv", ".tsv", ".parquet"],
        kind="area",
        label="Drag & drop a file",
    )
    file_browser = mo.ui.file_browser(
        filetypes=[".csv", ".tsv", ".parquet"],
        multiple=False,
        label="Or browse:",
    )
    cli_file_path = mo.cli_args().get("file", None)

    mo.vstack(
        [
            mo.md("# Capricho Data Explorer"),
            mo.md(
                "Load a Capricho aggregated output file (CSV, TSV, or Parquet) "
                "to interactively explore data quality, flag filtering, and binarization."
            ),
            mo.hstack([file_upload, file_browser], justify="start", gap=2),
        ]
    )
    return cli_file_path, file_browser, file_upload


@app.cell
def _(cli_file_path, file_browser, file_upload, io, load_dataframe, mo):
    _loaded_df = None
    _source_filename = None

    if cli_file_path is not None:
        _source_filename = str(cli_file_path)
        _loaded_df = load_dataframe(_source_filename, _source_filename)
    elif file_browser.value and len(file_browser.value) > 0:
        _path = file_browser.path()
        _source_filename = file_browser.name()
        _loaded_df = load_dataframe(_path, _source_filename)
    elif file_upload.value and len(file_upload.value) > 0:
        _source_filename = file_upload.name()
        _loaded_df = load_dataframe(io.BytesIO(file_upload.contents()), _source_filename)

    mo.stop(
        _loaded_df is None,
        mo.md("*Upload or browse to a Capricho data file to begin analysis.*"),
    )

    loaded_df = _loaded_df
    source_filename = _source_filename

    mo.md(f"Loaded **{source_filename}**")
    return loaded_df, source_filename


@app.cell
def _(detect_data_properties, loaded_df, mo, source_filename):
    data_props = detect_data_properties(loaded_df)

    _flags_str = ", ".join(data_props["present_flags"]) if data_props["present_flags"] else "none detected"
    _agg_str = "Yes (pipe-separated)" if data_props["is_aggregated"] else "No (single values per row)"

    mo.md(
        f"""
    ## Dataset Overview

    | Property | Value |
    |----------|-------|
    | **File** | {source_filename} |
    | **Rows** | {len(loaded_df):,} |
    | **Columns** | {len(loaded_df.columns)} |
    | **Aggregated** | {_agg_str} |
    | **Value column** | `{data_props["value_column"]}` |
    | **Axis limits** | {data_props["axis_limits"]} |
    | **Flags found** | {_flags_str} |
    """
    )
    return (data_props,)


@app.cell
def _(data_props, explode_assay_comparability, loaded_df, mo):
    mo.stop(
        not data_props["is_aggregated"],
        mo.md("*Comparability analysis requires aggregated data with pipe-separated values.*"),
    )

    _sep = "|"
    _val_col = data_props["value_column"]
    subset = loaded_df.query(f'{_val_col}.str.contains("{_sep}", regex=False)').assign(
        repeat=lambda x: range(len(x))
    )

    mo.stop(
        len(subset) == 0,
        mo.md("*No compounds with multiple assay measurements found.*"),
    )

    exploded_subset = explode_assay_comparability(subset, value_column=_val_col)

    mo.md(
        f"**{len(subset):,}** compounds with multiple measurements "
        f"→ **{len(exploded_subset):,}** pairwise assay comparisons"
    )
    return (exploded_subset,)


@app.cell
def _(
    data_props,
    exploded_subset,
    mo,
    np,
    plot_multi_panel_comparability,
    plt,
):
    mo.stop(not data_props["is_aggregated"])

    _val_col = data_props["value_column"]
    _flags = data_props["present_flags"]

    _ncols = min(len(_flags), 4) if _flags else 1
    _nrows = int(np.ceil(len(_flags) / _ncols)) if _flags else 1
    _figsize = (_ncols * 4.5, _nrows * 4.5)

    fig, axs = plot_multi_panel_comparability(
        exploded_subset,
        _flags,
        title="Comparability Across Data Quality Flags",
        figsize=_figsize,
        ncols=_ncols,
        value_column=_val_col,
        axis_limits=data_props["axis_limits"],
    )
    plt.gca()
    return


@app.cell
def _(data_props, mo):
    mo.stop(not data_props["is_aggregated"])

    mo.md("## Inspect Individual Flags")
    return


@app.cell
def _(data_props, mo):
    mo.stop(not data_props["is_aggregated"])

    _flags = data_props["present_flags"]
    mo.stop(len(_flags) == 0, mo.md("*No flags detected in this dataset.*"))

    comment_selector = mo.ui.dropdown(
        options=_flags,
        value=_flags[0],
        label="Data Quality Flag",
        full_width=True,
    )
    comment_selector
    return (comment_selector,)


@app.cell
def _(
    build_comparability_chart,
    build_query_string,
    comment_selector,
    data_props,
    exploded_subset,
    mo,
):
    _val_col = data_props["value_column"]
    _query_str = build_query_string(comment_selector.value, value_column=_val_col)
    filtered_data = exploded_subset.query(_query_str)

    mo.stop(
        len(filtered_data) == 0,
        mo.md(f"**No data available for:** {comment_selector.value}"),
    )

    filtered_data = filtered_data.copy()
    filtered_data[f"{_val_col}_x"] = filtered_data[f"{_val_col}_x"].astype(float)
    filtered_data[f"{_val_col}_y"] = filtered_data[f"{_val_col}_y"].astype(float)

    _chart_data = filtered_data[
        [f"{_val_col}_x", f"{_val_col}_y", "connectivity", "assay_chembl_id_x", "assay_chembl_id_y"]
    ].copy()

    _alt_chart = build_comparability_chart(
        _chart_data,
        title=f"Comparability: {comment_selector.value}",
        color="steelblue",
        value_column=_val_col,
        axis_limits=data_props["axis_limits"],
    )

    chart = mo.ui.altair_chart(_alt_chart)
    chart
    return chart, filtered_data


@app.cell
def _(chart, data_props, filtered_data, mo, smiles_to_svg_html):
    _val_col = data_props["value_column"]
    _selected_data = chart.apply_selection(filtered_data)

    if _selected_data is None or len(_selected_data) == 0:
        mo.md("*Click on points in the chart above to see details*")
        mo.stop(True)

    _data_to_render = _selected_data
    if len(_selected_data) > 50:
        mo.output.append(mo.md(f"**{len(_selected_data)} rows selected** (showing first 50 with structures)"))
        _data_to_render = _selected_data.head(50)

    _display_cols = [
        "connectivity",
        "canonical_smiles_x",
        "canonical_smiles_y",
        f"{_val_col}_x",
        f"{_val_col}_y",
        "assay_chembl_id_x",
        "assay_chembl_id_y",
        "data_processing_comment_x",
        "data_processing_comment_y",
        "data_dropping_comment_x",
        "data_dropping_comment_y",
    ]
    _available_cols = [c for c in _display_cols if c in _data_to_render.columns]
    _display_data = _data_to_render[_available_cols].copy()

    if "canonical_smiles_x" in _display_data.columns:
        _display_data["Molecule_x"] = _display_data["canonical_smiles_x"].apply(smiles_to_svg_html)
    if "canonical_smiles_y" in _display_data.columns:
        _display_data["Molecule_y"] = _display_data["canonical_smiles_y"].apply(smiles_to_svg_html)
    if "assay_chembl_id_x" in _display_data.columns:
        _display_data["Assay_x"] = _display_data["assay_chembl_id_x"].apply(
            lambda x: mo.Html(f'<a href="https://www.ebi.ac.uk/chembl/explore/assay/{x}" target="_blank">{x}</a>')
        )
    if "assay_chembl_id_y" in _display_data.columns:
        _display_data["Assay_y"] = _display_data["assay_chembl_id_y"].apply(
            lambda x: mo.Html(f'<a href="https://www.ebi.ac.uk/chembl/explore/assay/{x}" target="_blank">{x}</a>')
        )

    _table_cols = [
        "connectivity",
        "Molecule_x",
        "Molecule_y",
        f"{_val_col}_x",
        f"{_val_col}_y",
        "Assay_x",
        "Assay_y",
        "data_processing_comment_x",
        "data_processing_comment_y",
        "data_dropping_comment_x",
        "data_dropping_comment_y",
    ]
    _table_cols = [c for c in _table_cols if c in _display_data.columns]
    mo.ui.table(_display_data[_table_cols].to_dict("records"))
    return


@app.cell
def _(compute_correlation_metrics, data_props, filtered_data, mo):
    mo.stop(len(filtered_data) == 0)

    _val_col = data_props["value_column"]
    _metrics = compute_correlation_metrics(
        filtered_data[f"{_val_col}_x"],
        filtered_data[f"{_val_col}_y"],
    )

    mo.md(
        f"""
    ### Correlation Metrics
    - **R²**: {_metrics["r2"]:.3f}
    - **Spearman ρ**: {_metrics["spearman"]:.3f}
    - **Kendall τ**: {_metrics["kendall"]:.3f}
    """
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ## Data Preparation

    Select quality flags to exclude from the dataset. Flagged measurements are removed
    individually from aggregated rows (statistics are recalculated); rows are only dropped
    entirely if all their measurements are flagged.
    """)
    return


@app.cell
def _(DroppingComment, data_props, mo):
    _dropping_values = {dc.value for dc in DroppingComment}
    dropping_flag_list = [f for f in data_props["present_flags"] if any(f.startswith(dv) for dv in _dropping_values)]

    if len(dropping_flag_list) > 0:
        flag_checkboxes = mo.ui.array(
            [mo.ui.checkbox(label=f) for f in dropping_flag_list],
            label="Flags to exclude:",
        )
    else:
        flag_checkboxes = mo.ui.array([], label="No dropping flags found in data")

    deduplicate_checkbox = mo.ui.checkbox(label="Deduplicate identical values within aggregated rows")

    mo.vstack([flag_checkboxes, deduplicate_checkbox])
    return deduplicate_checkbox, dropping_flag_list, flag_checkboxes


@app.cell
def _(
    clean_data,
    data_props,
    deduplicate_checkbox,
    dropping_flag_list,
    flag_checkboxes,
    loaded_df,
    mo,
):
    _selected_flags = [
        dropping_flag_list[i] for i, checked in enumerate(flag_checkboxes.value) if checked
    ]
    _do_dedup = deduplicate_checkbox.value

    if len(_selected_flags) == 0 and not _do_dedup:
        prepared_df = loaded_df.copy()
        _output = mo.md("*Select flags to exclude or enable deduplication to preview cleaning.*")
    else:
        prepared_df = clean_data(
            loaded_df.copy(),
            drop_flags=_selected_flags if _selected_flags else None,
            deduplicate=_do_dedup,
            value_col=data_props["value_column"],
        )
        _pct = len(prepared_df) / len(loaded_df) * 100
        _excluded_str = ", ".join(_selected_flags) if _selected_flags else "none"
        _output = mo.md(
            f"""
    ### Preparation Result

    | | Count |
    |--|-------|
    | **Original rows** | {len(loaded_df):,} |
    | **After cleaning** | {len(prepared_df):,} |
    | **Retained** | {_pct:.1f}% |
    | **Flags excluded** | {_excluded_str} |
    | **Deduplicated** | {"Yes" if _do_dedup else "No"} |
    """
        )
    _output
    return (prepared_df,)


@app.cell
def _(
    build_comparability_chart,
    compute_correlation_metrics,
    data_props,
    deduplicate_checkbox,
    dropping_flag_list,
    explode_assay_comparability,
    exploded_subset,
    flag_checkboxes,
    mo,
    prepared_df,
):
    _selected_flags = [
        dropping_flag_list[i] for i, checked in enumerate(flag_checkboxes.value) if checked
    ]
    mo.stop(
        not data_props["is_aggregated"] or (len(_selected_flags) == 0 and not deduplicate_checkbox.value),
    )

    _val_col = data_props["value_column"]
    _sep = "|"

    _clean_subset = prepared_df.query(f'{_val_col}.str.contains("{_sep}", regex=False)').assign(
        repeat=lambda x: range(len(x))
    )

    mo.stop(len(_clean_subset) == 0, mo.md("*No multi-assay compounds remain after cleaning.*"))

    _clean_exploded = explode_assay_comparability(_clean_subset, value_column=_val_col)

    mo.stop(len(_clean_exploded) == 0)

    _clean_exploded[f"{_val_col}_x"] = _clean_exploded[f"{_val_col}_x"].astype(float)
    _clean_exploded[f"{_val_col}_y"] = _clean_exploded[f"{_val_col}_y"].astype(float)

    _orig_data = exploded_subset[
        [f"{_val_col}_x", f"{_val_col}_y", "connectivity", "assay_chembl_id_x", "assay_chembl_id_y"]
    ].copy()
    _orig_data[f"{_val_col}_x"] = _orig_data[f"{_val_col}_x"].astype(float)
    _orig_data[f"{_val_col}_y"] = _orig_data[f"{_val_col}_y"].astype(float)

    _clean_data = _clean_exploded[
        [f"{_val_col}_x", f"{_val_col}_y", "connectivity", "assay_chembl_id_x", "assay_chembl_id_y"]
    ].copy()

    _orig_chart = build_comparability_chart(
        _orig_data, "Before cleaning", "steelblue", _val_col, data_props["axis_limits"]
    )
    _clean_chart = build_comparability_chart(
        _clean_data, "After cleaning", "forestgreen", _val_col, data_props["axis_limits"]
    )

    _orig_metrics = compute_correlation_metrics(
        _orig_data[f"{_val_col}_x"], _orig_data[f"{_val_col}_y"]
    )
    _clean_metrics = compute_correlation_metrics(
        _clean_data[f"{_val_col}_x"], _clean_data[f"{_val_col}_y"]
    )

    mo.vstack(
        [
            mo.md("### Before vs After Cleaning"),
            mo.hstack(
                [mo.ui.altair_chart(_orig_chart), mo.ui.altair_chart(_clean_chart)],
                justify="center",
            ),
            mo.md(
                f"""
    | Metric | Before | After | Delta |
    |--------|--------|-------|-------|
    | **R²** | {_orig_metrics["r2"]:.3f} | {_clean_metrics["r2"]:.3f} | {_clean_metrics["r2"] - _orig_metrics["r2"]:+.3f} |
    | **Spearman ρ** | {_orig_metrics["spearman"]:.3f} | {_clean_metrics["spearman"]:.3f} | {_clean_metrics["spearman"] - _orig_metrics["spearman"]:+.3f} |
    | **Kendall τ** | {_orig_metrics["kendall"]:.3f} | {_clean_metrics["kendall"]:.3f} | {_clean_metrics["kendall"] - _orig_metrics["kendall"]:+.3f} |
    | **Pairs** | {len(_orig_data):,} | {len(_clean_data):,} | {len(_clean_data) - len(_orig_data):,} |
    """
            ),
        ]
    )
    return


@app.cell
def _(mo):
    mo.md("""
    ---
    ## Binarization

    Convert continuous activity values to binary labels (active/inactive) using a threshold.
    Adjust the threshold and conflict resolution strategy to preview the result.
    """)
    return


@app.cell
def _(VALID_CONFLICT_STRATEGIES, data_props, mo, pd, prepared_df):
    _val_col = data_props["value_column"]
    _mean_col = f"{_val_col}_mean"

    if _mean_col in prepared_df.columns:
        _vals = pd.to_numeric(prepared_df[_mean_col], errors="coerce").dropna()
        if len(_vals) > 0:
            _lo = max(float(_vals.min()) - 0.5, 0)
            _hi = float(_vals.max()) + 0.5
        else:
            _lo, _hi = 3.0, 12.0
    else:
        _lo, _hi = 3.0, 12.0

    _default = 6.0 if _val_col == "pchembl_value" else (_lo + _hi) / 2

    threshold_slider = mo.ui.slider(
        start=_lo,
        stop=_hi,
        value=_default,
        step=0.1,
        label="Activity threshold",
        full_width=True,
    )

    conflict_dropdown = mo.ui.dropdown(
        options=["None"] + sorted(VALID_CONFLICT_STRATEGIES),
        value="None",
        label="Conflict resolution",
    )

    mo.hstack([threshold_slider, conflict_dropdown], justify="start", gap=2)
    return conflict_dropdown, threshold_slider


@app.cell
def _(
    alt,
    binarize_aggregated_data,
    conflict_dropdown,
    data_props,
    mo,
    pd,
    prepared_df,
    threshold_slider,
):
    _val_col = data_props["value_column"]
    _mean_col = f"{_val_col}_mean"

    mo.stop(
        _mean_col not in prepared_df.columns,
        mo.md(f"*Binarization requires a `{_mean_col}` column. Run `capricho get` to produce aggregated data.*"),
    )

    _conflict_res = None if conflict_dropdown.value == "None" else conflict_dropdown.value

    binarized_df = binarize_aggregated_data(
        prepared_df.copy(),
        threshold=threshold_slider.value,
        value_column=_mean_col,
        conflict_resolution=_conflict_res,
    )

    _active = (binarized_df["activity_binary"] == 1).sum()
    _inactive = (binarized_df["activity_binary"] == 0).sum()
    _unclassified = binarized_df["activity_binary"].isna().sum()
    _total = len(binarized_df)

    _dist_data = pd.DataFrame(
        {
            "Label": ["Active (1)", "Inactive (0)", "Unclassified (NaN)"],
            "Count": [_active, _inactive, _unclassified],
        }
    )

    _bar_chart = (
        alt.Chart(_dist_data)
        .mark_bar()
        .encode(
            x=alt.X("Label:N", title=None, sort=["Active (1)", "Inactive (0)", "Unclassified (NaN)"]),
            y=alt.Y("Count:Q", title="Count"),
            color=alt.Color(
                "Label:N",
                scale=alt.Scale(
                    domain=["Active (1)", "Inactive (0)", "Unclassified (NaN)"],
                    range=["#2ca02c", "#d62728", "#7f7f7f"],
                ),
                legend=None,
            ),
            tooltip=["Label:N", "Count:Q"],
        )
        .properties(width=400, height=300, title=f"Binarization at threshold = {threshold_slider.value:.1f}")
    )

    mo.vstack(
        [
            mo.hstack([mo.ui.altair_chart(_bar_chart)], justify="center"),
            mo.md(
                f"""
    ### Binarization Summary (threshold = {threshold_slider.value:.1f})

    | | Count | % |
    |--|-------|---|
    | **Active** | {_active:,} | {_active / _total * 100:.1f}% |
    | **Inactive** | {_inactive:,} | {_inactive / _total * 100:.1f}% |
    | **Unclassified** | {_unclassified:,} | {_unclassified / _total * 100:.1f}% |
    | **Total** | {_total:,} | |
    """
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
