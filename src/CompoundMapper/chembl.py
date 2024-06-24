"""Module holding functionalities for the ChEMBL API."""

from typing import List, Tuple, Union

import numpy as np

import pandas as pd
from chembl_webresource_client.new_client import new_client
from .logger import logger

assays_api = new_client.assay
activity_api = new_client.activity
compounds_api = new_client.molecule

# temporary; to be removed after pandas 3.0
# check for pandas' version. If lower than 3.0, set the option to avoid silent downcasting
if pd.__version__ < "3.0.0":
    pd.set_option("future.no_silent_downcasting", True)

# Info on Chirality:
# The chirality flag shows whether a drug is dosed as a racemic mixture (0), single stereoisomer (1) or as an achiral molecule (2), for unchecked compounds the chirality flag = -1.
# source: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/drug-and-compound-questions#:~:text=Blog%20post.-,Can%20you%20provide%20more%20details%20on%20the%20chirality%20flag%3F,-The%20chirality%20flag


def find_dict_in_dataframe(df):
    cols_w_dicts = []
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, dict)).any():
            logger.info(f"Column '{col}' contains dictionaries.")
            logger.info(
                "Rows with dictionaries: "
                " ".join(
                    df[df[col].apply(lambda x: isinstance(x, dict))].index.astype(str)
                )
            )
            cols_w_dicts.append(col)
    if cols_w_dicts:
        return cols_w_dicts


def get_publications_details(assay_chembl_ids: list) -> dict:
    """From a list of ChEMBL assay IDs, get the publication details.
    Args:
        assay_chembl_ids: list of ChEMBL assay IDs.
    Returns:
        dict: a dictionary with assay IDs as keys and publication details as values.
    """
    assay = new_client.assay
    document = new_client.document

    publications_details = {}
    for assay_id in assay_chembl_ids:
        assay_data = assay.get(assay_id)
        authors, doi, journal, volume, year, title = None, None, None, None, None, None
        chembl_release = None
        if assay_data and "document_chembl_id" in assay_data:
            doc_data = document.get(assay_data["document_chembl_id"])
            if doc_data:
                authors = doc_data.get("authors")
                doi = doc_data.get("doi")
                journal = doc_data.get("journal")
                volume = doc_data.get("volume")
                year = doc_data.get("year")
                title = doc_data.get("title")
                chembl_release = doc_data.get("chembl_release")
                if isinstance(chembl_release, dict):
                    chembl_release = chembl_release.get("chembl_release")

        publications_details[assay_id] = {
            "Authors": authors,
            "DOI": doi,
            "Journal": journal,
            "Volume": volume,
            "Year": year,
            "Title": title,
            "ChEMBL_Release": chembl_release,
        }
    return publications_details


def molecule_info_from_chembl(molecule_chembl_id: list) -> dict:
    """Get information on a molecule from ChEMBL.
    Args:
        molecule_chembl_id: a list of molecule ChEMBL IDs.
    Returns:
        pd.DataFrame: a DataFrame with the molecule information.
    """
    extracted = {}
    result = compounds_api.filter(molecule_chembl_id__in=molecule_chembl_id).only(
        "molecule_hierarchy",
        "molecule_structures",
        "chemical_probes",
        "chirality",
        "oral",
        "prodrug",
        "max_phase",
        "therapeutical_flag",
        "withdrawn_flag",
        "indication_class",
    )
    if result:
        for r, mol_id in zip(result, molecule_chembl_id):
            if r is None:
                logger.warning(f"No information found for molecule {mol_id}")
                continue
            if r["molecule_hierarchy"] is not None:
                hierarchy_active_id = r.get("molecule_hierarchy", {}).get(
                    "active_chembl_id", None
                )
                hierarchy_molecule_id = r.get("molecule_hierarchy", {}).get(
                    "molecule_chembl_id", None
                )
                hierarchy_parent_id = r.get("molecule_hierarchy", {}).get(
                    "parent_chembl_id", None
                )
                r.pop("molecule_hierarchy")
            else:
                logger.warning(f"No hierarchy information found for molecule {mol_id}")
                hierarchy_active_id = None
                hierarchy_molecule_id = None
                hierarchy_parent_id = None
            if r["molecule_structures"] is not None:
                canonical_smiles = r.get("molecule_structures", {}).get(
                    "canonical_smiles", None
                )
                standard_inchikey = r.get("molecule_structures", {}).get(
                    "standard_inchi_key", None
                )
                r.pop("molecule_structures")
            else:
                logger.warning(f"No structure information found for molecule {mol_id}")
                canonical_smiles = None
                standard_inchikey = None
            extracted[mol_id] = {
                "hierarchy_active_id": hierarchy_active_id,
                "hierarchy_molecule_id": hierarchy_molecule_id,
                "hierarchy_parent_id": hierarchy_parent_id,
                "canonical_smiles": canonical_smiles,
                "standard_inchikey": standard_inchikey,
                **r,
            }
    else:
        raise ValueError(f"No information found for {molecule_chembl_id}")
    return pd.DataFrame.from_dict(extracted, orient="index").reset_index(
        names="molecule_chembl_id"
    )


def assay_info_from_chembl(
    assay_chembl_ids: list,
    confidence_scores: list | None = None,
    **kwargs,
) -> pd.DataFrame:
    """Take a list of assay chembl ids and get their respective assays in ChEMBL.
    Args:
        assay_chembl_ids: a list of assay ChEMBL IDs.
        kwargs: keywords arguments to filter the assays.
    Returns:
        pd.DataFrame: a DataFrame with the assays.
    """
    if confidence_scores is None:
        confidence_scores = list(range(0, 10))
    activity_kwargs = {
        "assay_chembl_id__in": assay_chembl_ids,
        "confidence_score__in": confidence_scores,
        **kwargs,
    }
    assays = assays_api.filter(**activity_kwargs).only(
        "assay_cell_type",
        "assay_chembl_id",
        "assay_organism",
        "assay_subcellular_fraction",
        "assay_tissue",
        "assay_type",
        "assay_type_description",
        "confidence_score",
        "description",
        "document_chembl_id",
        "relationship_description",
        "relationship_type",
        "target_chembl_id",
        "variant_sequence",
    )
    if assays:
        assays_df = pd.DataFrame.from_records(assays)
        if find_dict_in_dataframe(assays_df) is not None:
            logger.warning("Keeping only mutation info from `variant_sequence`.")
            assays_df = assays_df.assign(
                variant_sequence=lambda x: x.variant_sequence.apply(
                    lambda y: y.get("mutation") if isinstance(y, dict) else y
                )
            )
    else:
        activity_kwargs.pop("assay_chembl_id__in")
        raise ValueError(
            f"No assays found for the ids: {assay_chembl_ids} with the parameters: {activity_kwargs}"
        )
    return assays_df


def bioactivities_from_chembl(
    molecule_chembl_ids: list,
    **kwargs,
) -> pd.DataFrame:
    """Take a list of molecule chembl ids and get their respective bioactivities in ChEMBL.
    Args:
        molecule_chembl_id: list of molecule ChEMBL IDs.
        kwargs: example -> `standard_relation="=", assay_type__in=["B", "F"]`.
    Returns:
        pd.DataFrame: a DataFrame with the bioactivities.
    """
    activity_kwargs = {
        "molecule_chembl_id__in": molecule_chembl_ids,
        **kwargs,
    }
    bioactivities = activity_api.filter(**activity_kwargs).only(
        "activity_id",
        "assay_chembl_id",
        "assay_description",
        "assay_type",
        "molecule_chembl_id",
        "standard_flag",
        "standard_relation",
        "standard_type",
        "standard_units",
        "standard_value",
        "pchembl_value",
        "target_chembl_id",
        "target_organism",
        "data_validity_description",
        "potential_duplicate",
    )
    if bioactivities:
        assays_df = pd.DataFrame.from_records(bioactivities)
    else:
        raise ValueError(f"No bioactivities found for the ids: {molecule_chembl_ids}")
    return assays_df


def convert_to_log10(df):
    """Function to be applied to the whole DataFrame. Will convert the standard_value
    column to pchembl_value column, if the standard_units are in nM, µM or uM.
    Args:
        df: a bioactivity DataFrame. e.g.: output from `bioactivities_from_chembl`.
    """

    def compute_log(row):
        value = row["standard_value"]
        unit = row["standard_units"]
        if unit == "µM" or unit == "uM":
            value_in_M = value * 1e-6
        elif unit == "nM":
            value_in_M = value * 1e-9
        return -np.log10(value_in_M)

    desired_units = ["nM", "µM", "uM"]
    final_df = []
    for unit in desired_units:  # noqa: B007
        tmp_df = df.query("standard_units == @unit")
        if tmp_df.shape[0] > 0:
            tmp_df = tmp_df.assign(pchembl_value=lambda x: x.apply(compute_log, axis=1))
            final_df.append(tmp_df)
        else:
            continue
    if final_df == []:
        return None
    else:
        return pd.concat(final_df, ignore_index=True)


def process_bioactivities(bioactivities_df: pd.DataFrame) -> pd.DataFrame:
    """Processes the bioactivities DataFrame. Will convert the standard_value
    column to pchembl_value column if the standard_units are in nM, µM or uM.
    Args:
        bioactivities_df: bioactivity dataframe, e.g.: output from `bioactivities_from_chembl`.
    Returns:
        pd.DataFrame: the processed bioactivities DataFrame.
    """
    bioactivities_df = (
        bioactivities_df.astype(
            {"standard_value": "float32", "pchembl_value": "float32"}
        )
        .replace({None: np.nan})
        .query("data_validity_description.isna()")
        .query("potential_duplicate == 0")
        .drop(columns=["data_validity_description", "potential_duplicate"])
    )
    with_pchembl = bioactivities_df.query("~pchembl_value.isna()")
    without_pchembl = convert_to_log10(  # if pchembl value not present, calculate it
        bioactivities_df.query("pchembl_value.isna()")
    )
    if without_pchembl is not None:
        bioactivities_df = pd.concat([with_pchembl, without_pchembl])
    else:
        bioactivities_df = with_pchembl
    return bioactivities_df


def filtered_data_workflow(
    molecule_ids: list, confidence_scores: Union[List, Tuple] = (9, 8)
) -> pd.DataFrame:
    """
    Retrieves and merges data from ChEMBL for given molecule IDs, filtered by specified confidence scores.
    Args:
        molecule_ids: List of molecule ChEMBL IDs.
        confidence_scores: List of confidence scores to filter the data. Defaults to [9, 8].
    Returns:
        pd.DataFrame: Merged DataFrame with molecule, bioactivity, and assay information.
    """

    # Step 1: Get bioactivities for given molecule IDs
    bioactivities_df = process_bioactivities(
        bioactivities_from_chembl(
            molecule_ids, standard_relation="=", assay_type__in=["B", "F"]
        )
    )
    # Step 2: Extract unique assay IDs from the bioactivities DataFrame
    unique_assay_ids = bioactivities_df["assay_chembl_id"].unique().tolist()
    # Step 3: Get assay information for these unique assay IDs & merge data
    assays_df = assay_info_from_chembl(
        unique_assay_ids, confidence_scores=list(confidence_scores)
    )
    merged_df = pd.merge(
        bioactivities_df,
        assays_df,
        on=["assay_chembl_id", "target_chembl_id", "assay_type"],
        how="inner",
    ).drop_duplicates()
    return merged_df
