import marimo

__generated_with = "0.17.2"
app = marimo.App(width="medium")


@app.cell
def _():
    import altair as alt
    import marimo as mo
    import matplotlib.pyplot as plt
    import pandas as pd
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D

    from Capricho.analysis import (
        build_query_string,
        explode_assay_comparability,
        get_all_comments,
        plot_multi_panel_comparability,
    )

    return (
        AllChem,
        Chem,
        alt,
        build_query_string,
        explode_assay_comparability,
        get_all_comments,
        mo,
        pd,
        plot_multi_panel_comparability,
        plt,
        rdMolDraw2D,
    )


@app.cell
def _(AllChem, Chem, mo, pd, rdMolDraw2D):
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

    return (smiles_to_svg_html,)


@app.cell
def _(mo):
    mo.md("""# IC50 Data Quality Analysis""")
    return


@app.cell
def _(mo):
    mo.md(
        """
    This interactive notebook allows you to explore how different data quality flags
    affect the comparability of pChEMBL measurements across assays. Click on points
    in the scatter plot to see the underlying data.
    """
    )
    return


@app.cell
def _(mo, pd):
    df = pd.read_csv("curated_data/curated-IC50-NoAssayOverlap.csv", engine="pyarrow")
    subset = df.query('processed_smiles.str.contains("|", regex=False)').assign(
        repeat=lambda x: range(len(x))
    )
    mo.md(f"Loaded **{len(df):,}** total compounds, **{len(subset):,}** have multiple assay measurements")
    return (subset,)


@app.cell
def _(explode_assay_comparability, get_all_comments, mo, subset):
    exploded_subset = explode_assay_comparability(subset)
    all_comments = get_all_comments()
    mo.md(f"Created **{len(exploded_subset):,}** pairwise assay comparisons for analysis")
    return all_comments, exploded_subset


@app.cell
def _(mo):
    mo.md("""## Overview: Comparability Across All Data Quality Flags""")
    return


@app.cell
def _(all_comments, exploded_subset, plot_multi_panel_comparability, plt):
    fig, axs = plot_multi_panel_comparability(
        exploded_subset,
        all_comments,
        title="IC50 Data Comparability by Quality Flag",
        figsize=(20, 8),
        ncols=5,
    )
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md("""## Select a data quality flag to analyze:""")
    return


@app.cell
def _(all_comments, mo):
    comment_selector = mo.ui.dropdown(
        options=all_comments,
        value=all_comments[0],
        label="Data Quality Flag",
        full_width=True,
    )
    comment_selector
    return (comment_selector,)


@app.cell
def _(build_query_string, comment_selector, exploded_subset, mo):
    query_str = build_query_string(comment_selector.value)
    filtered_data = exploded_subset.query(query_str)

    if len(filtered_data) == 0:
        mo.md(f"**No data available for:** {comment_selector.value}")
    else:
        filtered_data = filtered_data.copy()
        filtered_data["pchembl_value_x"] = filtered_data["pchembl_value_x"].astype(float)
        filtered_data["pchembl_value_y"] = filtered_data["pchembl_value_y"].astype(float)
        mo.md(f"**{len(filtered_data):,}** pairwise comparisons with this flag")
    return (filtered_data,)


@app.cell
def _(alt, comment_selector, filtered_data, mo, pd):
    mo.stop(len(filtered_data) == 0)

    chart_data = filtered_data[
        ["pchembl_value_x", "pchembl_value_y", "connectivity", "assay_chembl_id_x", "assay_chembl_id_y"]
    ].copy()

    base_chart = (
        alt.Chart(chart_data)
        .mark_circle(size=60, opacity=0.5, color="steelblue")
        .encode(
            x=alt.X("pchembl_value_x:Q", title="Assay 1 pChEMBL value", scale=alt.Scale(domain=[3, 12])),
            y=alt.Y("pchembl_value_y:Q", title="Assay 2 pChEMBL value", scale=alt.Scale(domain=[3, 12])),
            tooltip=[
                alt.Tooltip("connectivity:N", title="Compound"),
                alt.Tooltip("pchembl_value_x:Q", title="Assay 1 pChEMBL", format=".2f"),
                alt.Tooltip("pchembl_value_y:Q", title="Assay 2 pChEMBL", format=".2f"),
                alt.Tooltip("assay_chembl_id_x:N", title="Assay 1 ID"),
                alt.Tooltip("assay_chembl_id_y:N", title="Assay 2 ID"),
            ],
        )
    )

    identity_line = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [3, 12]}))
        .mark_line(color="black", strokeDash=[0])
        .encode(x="x:Q", y="y:Q")
    )

    plus_one = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [4, 13]}))
        .mark_line(color="black", strokeDash=[5, 5])
        .encode(x="x:Q", y="y:Q")
    )

    minus_one = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [2, 11]}))
        .mark_line(color="black", strokeDash=[5, 5])
        .encode(x="x:Q", y="y:Q")
    )

    plus_03 = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [3.3, 12.3]}))
        .mark_line(color="black", strokeDash=[2, 2])
        .encode(x="x:Q", y="y:Q")
    )

    minus_03 = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [2.7, 11.7]}))
        .mark_line(color="black", strokeDash=[2, 2])
        .encode(x="x:Q", y="y:Q")
    )

    _chart = (base_chart + identity_line + plus_one + minus_one + plus_03 + minus_03).properties(
        width=600, height=600, title=f"IC50 Comparability: {comment_selector.value}"
    )

    chart = mo.ui.altair_chart(_chart)
    chart
    return (chart,)


@app.cell
def _(mo):
    mo.md("""### Selected Data Points""")
    return


@app.cell
def _(chart, filtered_data, mo, smiles_to_svg_html):
    _selected_data = chart.apply_selection(filtered_data)

    if _selected_data is None or len(_selected_data) == 0:
        mo.md("*Click on points in the chart above to see details*")
        mo.stop(True)

    _data_to_render = _selected_data
    if len(_selected_data) > 50:
        mo.output.append(mo.md(f"**{len(_selected_data)} rows selected** (showing first 50 with structures)"))
        _data_to_render = _selected_data.head(50)

    _display_data = _data_to_render[
        [
            "connectivity",
            "canonical_smiles_x",
            "canonical_smiles_y",
            "pchembl_value_x",
            "pchembl_value_y",
            "assay_chembl_id_x",
            "assay_chembl_id_y",
            "data_processing_comment_x",
            "data_processing_comment_y",
            "data_dropping_comment_x",
            "data_dropping_comment_y",
        ]
    ].copy()

    _display_data["Molecule_x"] = _display_data["canonical_smiles_x"].apply(smiles_to_svg_html)
    _display_data["Molecule_y"] = _display_data["canonical_smiles_y"].apply(smiles_to_svg_html)
    _display_data["Assay_x"] = _display_data["assay_chembl_id_x"].apply(
        lambda x: mo.Html(f'<a href="https://www.ebi.ac.uk/chembl/explore/assay/{x}" target="_blank">{x}</a>')
    )
    _display_data["Assay_y"] = _display_data["assay_chembl_id_y"].apply(
        lambda x: mo.Html(f'<a href="https://www.ebi.ac.uk/chembl/explore/assay/{x}" target="_blank">{x}</a>')
    )

    _table_data = _display_data[
        [
            "connectivity",
            "Molecule_x",
            "Molecule_y",
            "pchembl_value_x",
            "pchembl_value_y",
            "Assay_x",
            "Assay_y",
            "data_processing_comment_x",
            "data_processing_comment_y",
            "data_dropping_comment_x",
            "data_dropping_comment_y",
        ]
    ].to_dict("records")
    mo.ui.table(_table_data)
    return


@app.cell
def _(filtered_data, mo):
    mo.stop(len(filtered_data) == 0)

    from scipy import stats
    from sklearn.metrics import r2_score

    xp = filtered_data["pchembl_value_x"]
    yp = filtered_data["pchembl_value_y"]

    r, _ = stats.spearmanr(xp, yp)
    tau, _ = stats.kendalltau(xp, yp)
    r2 = r2_score(xp, yp)

    mo.md(
        f"""
    ### Correlation Metrics
    - **R²**: {r2:.3f}
    - **Spearman ρ**: {r:.3f}
    - **Kendall τ**: {tau:.3f}
    """
    )
    return r2_score, stats


@app.cell
def _(mo):
    mo.md(
        """
    ---
    ## Clean Data Analysis
    Select which quality flags you want to exclude from the analysis.
    The plot below will show assay comparability after removing compounds with the selected flags.
    """
    )
    return


@app.cell
def _(all_comments, mo):
    flag_checkboxes = mo.ui.array(
        [mo.ui.checkbox(label=comment) for comment in all_comments], label="Flags to exclude:"
    )
    flag_checkboxes
    return (flag_checkboxes,)


@app.cell
def _(all_comments, exploded_subset, flag_checkboxes, mo, pd):
    selected_flags = [all_comments[i] for i in range(len(flag_checkboxes.value)) if flag_checkboxes.value[i]]

    if len(selected_flags) == 0:
        mo.md("*Select at least one flag to exclude to see the filtered data*")
        clean_data = pd.DataFrame()
    else:
        clean_data = exploded_subset.copy()

        for flag in selected_flags:
            mask = (
                ~clean_data["data_processing_comment_x"].fillna("-//-").str.contains(flag, regex=False)
                & ~clean_data["data_processing_comment_y"].fillna("-//-").str.contains(flag, regex=False)
                & ~clean_data["data_dropping_comment_x"].fillna("-//-").str.contains(flag, regex=False)
                & ~clean_data["data_dropping_comment_y"].fillna("-//-").str.contains(flag, regex=False)
            )
            clean_data = clean_data[mask]

        if len(clean_data) == 0:
            mo.md(f"**No data remaining after excluding:** {', '.join(selected_flags)}")
        else:
            clean_data = clean_data.copy()
            clean_data["pchembl_value_x"] = clean_data["pchembl_value_x"].astype(float)
            clean_data["pchembl_value_y"] = clean_data["pchembl_value_y"].astype(float)

            mo.md(
                f"""
        **Excluded flags:** {', '.join(selected_flags)}

        **{len(clean_data):,}** pairwise comparisons remaining ({len(clean_data)/len(exploded_subset)*100:.1f}% of original data)
        """
            )
    return (clean_data,)


@app.cell
def _(alt, clean_data, mo, pd):
    mo.stop(len(clean_data) == 0)

    clean_chart_data = clean_data[
        ["pchembl_value_x", "pchembl_value_y", "connectivity", "assay_chembl_id_x", "assay_chembl_id_y"]
    ].copy()

    clean_base = (
        alt.Chart(clean_chart_data)
        .mark_circle(size=60, opacity=0.5, color="forestgreen")
        .encode(
            x=alt.X("pchembl_value_x:Q", title="Assay 1 pChEMBL value", scale=alt.Scale(domain=[3, 12])),
            y=alt.Y("pchembl_value_y:Q", title="Assay 2 pChEMBL value", scale=alt.Scale(domain=[3, 12])),
            tooltip=[
                alt.Tooltip("connectivity:N", title="Compound"),
                alt.Tooltip("pchembl_value_x:Q", title="Assay 1 pChEMBL", format=".2f"),
                alt.Tooltip("pchembl_value_y:Q", title="Assay 2 pChEMBL", format=".2f"),
                alt.Tooltip("assay_chembl_id_x:N", title="Assay 1 ID"),
                alt.Tooltip("assay_chembl_id_y:N", title="Assay 2 ID"),
            ],
        )
    )

    clean_identity = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [3, 12]}))
        .mark_line(color="black", strokeDash=[0])
        .encode(x="x:Q", y="y:Q")
    )

    clean_plus_one = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [4, 13]}))
        .mark_line(color="black", strokeDash=[5, 5])
        .encode(x="x:Q", y="y:Q")
    )

    clean_minus_one = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [2, 11]}))
        .mark_line(color="black", strokeDash=[5, 5])
        .encode(x="x:Q", y="y:Q")
    )

    clean_plus_03 = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [3.3, 12.3]}))
        .mark_line(color="black", strokeDash=[2, 2])
        .encode(x="x:Q", y="y:Q")
    )

    clean_minus_03 = (
        alt.Chart(pd.DataFrame({"x": [3, 12], "y": [2.7, 11.7]}))
        .mark_line(color="black", strokeDash=[2, 2])
        .encode(x="x:Q", y="y:Q")
    )

    _clean_chart = (
        clean_base + clean_identity + clean_plus_one + clean_minus_one + clean_plus_03 + clean_minus_03
    ).properties(width=600, height=600, title="IC50 Comparability: Clean Data")

    clean_chart = mo.ui.altair_chart(_clean_chart)
    clean_chart
    return (clean_chart,)


@app.cell
def _(mo):
    mo.md("""### Selected Clean Data Points""")
    return


@app.cell
def _(clean_chart, clean_data, mo, smiles_to_svg_html):
    _selected_clean = clean_chart.apply_selection(clean_data)

    if _selected_clean is None or len(_selected_clean) == 0:
        mo.md("*Click on points in the chart above to see details*")
        mo.stop(True)

    _data_to_render = _selected_clean
    if len(_selected_clean) > 50:
        mo.output.append(
            mo.md(f"**{len(_selected_clean)} rows selected** (showing first 50 with structures)")
        )
        _data_to_render = _selected_clean.head(50)

    _display_clean = _data_to_render[
        [
            "connectivity",
            "canonical_smiles_x",
            "canonical_smiles_y",
            "pchembl_value_x",
            "pchembl_value_y",
            "assay_chembl_id_x",
            "assay_chembl_id_y",
            "data_processing_comment_x",
            "data_processing_comment_y",
            "data_dropping_comment_x",
            "data_dropping_comment_y",
        ]
    ].copy()

    _display_clean["Molecule_x"] = _display_clean["canonical_smiles_x"].apply(smiles_to_svg_html)
    _display_clean["Molecule_y"] = _display_clean["canonical_smiles_y"].apply(smiles_to_svg_html)
    _display_clean["Assay_x"] = _display_clean["assay_chembl_id_x"].apply(
        lambda x: mo.Html(f'<a href="https://www.ebi.ac.uk/chembl/explore/assay/{x}" target="_blank">{x}</a>')
    )
    _display_clean["Assay_y"] = _display_clean["assay_chembl_id_y"].apply(
        lambda x: mo.Html(f'<a href="https://www.ebi.ac.uk/chembl/explore/assay/{x}" target="_blank">{x}</a>')
    )

    _table_data_clean = _display_clean[
        [
            "connectivity",
            "Molecule_x",
            "Molecule_y",
            "pchembl_value_x",
            "pchembl_value_y",
            "Assay_x",
            "Assay_y",
            "data_processing_comment_x",
            "data_processing_comment_y",
            "data_dropping_comment_x",
            "data_dropping_comment_y",
        ]
    ].to_dict("records")
    mo.ui.table(_table_data_clean)
    return


@app.cell
def _(clean_data, mo, r2_score, stats):
    mo.stop(len(clean_data) == 0)

    xp_clean = clean_data["pchembl_value_x"]
    yp_clean = clean_data["pchembl_value_y"]

    r_clean, _ = stats.spearmanr(xp_clean, yp_clean)
    tau_clean, _ = stats.kendalltau(xp_clean, yp_clean)
    r2_clean = r2_score(xp_clean, yp_clean)

    mo.md(
        f"""
    ### Clean Data Correlation Metrics
    - **R²**: {r2_clean:.3f}
    - **Spearman ρ**: {r_clean:.3f}
    - **Kendall τ**: {tau_clean:.3f}
    """
    )
    return


if __name__ == "__main__":
    app.run()
