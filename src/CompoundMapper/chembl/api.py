"""Module holding functionalities for the ChEMBL API."""

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional

import pandas as pd
from chembl_webresource_client.new_client import new_client

from ..core.pandas_helper import find_dict_in_dataframe
from ..logger import logger
from .rate_limit import rate_limit

# Info on Chirality:
# The chirality flag shows whether a drug is dosed as a racemic mixture (0), single stereoisomer (1) or as an achiral molecule (2), for unchecked compounds the chirality flag = -1.
# source: https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/drug-and-compound-questions#:~:text=Blog%20post.-,Can%20you%20provide%20more%20details%20on%20the%20chirality%20flag%3F,-The%20chirality%20flag


def get_document_details(document_chembl_ids: list, chembl_version: Optional[int] = None) -> pd.DataFrame:
    """From a list of ChEMBL assay IDs, get the publication details.
    Args:
        assay_chembl_ids: list of ChEMBL assay IDs.
    Returns:
        dict: a dictionary with assay IDs as keys and publication details as values.
    """
    query_kwargs = {}
    if document_chembl_ids is not None:
        query_kwargs.update({"document_chembl_ids__in": document_chembl_ids})

    document_api = new_client.document
    documents = document_api.filter(**query_kwargs).only(
        "document_chembl_id",
        "doc_type",
        "authors",
        "doi",
        "journal",
        "volume",
        "year",
        "title",
        "chembl_release",
    )
    publications_details = {}
    for doc_data in documents:
        authors, doi, journal, volume, year, title = None, None, None, None, None, None
        chembl_release = None
        if doc_data:
            document_id = doc_data.get("document_chembl_id")
            doc_type = doc_data.get("doc_type")
            authors = doc_data.get("authors")
            doi = doc_data.get("doi")
            journal = doc_data.get("journal")
            volume = doc_data.get("volume")
            year = doc_data.get("year")
            title = doc_data.get("title")
            chembl_release = doc_data.get("chembl_release")
            if isinstance(chembl_release, dict):
                chembl_release = chembl_release.get("chembl_release")

            if chembl_version is not None:
                if chembl_release is None or int(chembl_release.split("_")[1]) > chembl_version:
                    continue

        publications_details[document_id] = {
            "doc_type": doc_type,
            "authors": authors,
            "doi": doi,
            "journal": journal,
            "volume": volume,
            "year": year,
            "title": title,
            "chembl_release": chembl_release,
        }
    return pd.DataFrame.from_dict(publications_details, orient="index").reset_index(
        names=["document_chembl_id"]
    )


def molecule_info_from_chembl(molecule_chembl_ids: list) -> pd.DataFrame:
    """Get information on a molecule from ChEMBL.
    Args:
        molecule_chembl_ids: a list of molecule ChEMBL IDs.
    Returns:
        pd.DataFrame: a DataFrame with the molecule information.
    """
    extracted = {}
    compounds_api = new_client.molecule
    result = compounds_api.filter(molecule_chembl_id__in=molecule_chembl_ids).only(
        "molecule_chembl_id",
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
        for idx, r in enumerate(result):
            mol_id = r.get("molecule_chembl_id")
            if r is None:
                logger.warning(f"No information found for molecule {mol_id}")
                continue
            extracted[idx] = parse_molecule_response(r, mol_id)
    else:
        raise ValueError(f"No information found for {molecule_chembl_ids}")
    return (
        pd.DataFrame.from_dict(extracted, orient="index")
        .sort_values(by="molecule_chembl_id", key=lambda col: col.map(lambda e: molecule_chembl_ids.index(e)))
        .reset_index(drop=True)
    )


def parse_molecule_response(r: dict, compound_id: str) -> dict:
    """Parse the response from the ChEMBL API for a molecule.

    Args:
        r: response, a dictionary with the information on the compound.
        compound_id: identifier to log warnings for the compound when no information is found.

    Returns:
        dict: a dictionary with the parsed information.
    """
    if r["molecule_hierarchy"] is not None:
        hierarchy_active_id = r.get("molecule_hierarchy", {}).get("active_chembl_id", None)
        hierarchy_molecule_id = r.get("molecule_hierarchy", {}).get("molecule_chembl_id", None)
        hierarchy_parent_id = r.get("molecule_hierarchy", {}).get("parent_chembl_id", None)
        r.pop("molecule_hierarchy")
    else:
        logger.warning(f"No hierarchy information found for compound {compound_id}")
        hierarchy_active_id = None
        hierarchy_molecule_id = None
        hierarchy_parent_id = None
    if r["molecule_structures"] is not None:
        canonical_smiles = r.get("molecule_structures", {}).get("canonical_smiles", None)
        standard_inchikey = r.get("molecule_structures", {}).get("standard_inchi_key", None)
        r.pop("molecule_structures")
    else:
        logger.warning(f"No structure information found for compound {compound_id}")
        canonical_smiles = None
        standard_inchikey = None
    return {
        "hierarchy_active_id": hierarchy_active_id,
        "hierarchy_molecule_id": hierarchy_molecule_id,
        "hierarchy_parent_id": hierarchy_parent_id,
        "canonical_smiles": canonical_smiles,
        "standard_inchikey": standard_inchikey,
        **r,
    }


def get_similar_compounds(smi: str, similarity: float) -> pd.DataFrame:
    """Fetch similar compounds from ChEMBL using the similarity API.

    Args:
        smiles: single smiles string to find similar molecules to.
        similarity: similarity threshold to use for the search.

    Returns:
        pd.DataFrame: a DataFrame with the similar molecules.
    """
    similarity_api = new_client.similarity
    extracted = {}
    result = similarity_api.filter(smiles=smi, similarity=similarity).only(
        "molecule_chembl_id",
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
        "similarity",
    )
    if result:
        for idx, r in enumerate(result):
            if r is None:
                logger.warning(f"No similar molecules were found for reponse {idx}")
            else:
                extracted[idx] = {"querySmiles": smi, **parse_molecule_response(r, smi)}
    else:
        logger.warning(f"No similar molecules were found to {smi}")
    return pd.DataFrame.from_dict(extracted, orient="index")


@rate_limit(max_per_second=5)
def get_similars_from_smiles(smiles: list[str], similarity: float, n_threads: int = 1) -> pd.DataFrame:
    """Use the ChEMBL API to get similar compounds to a list of SMILES. Though multiple threads
    can be used, the rate limit is set to 5 calls per second not to overload the API.

    Args:
        smiles: list of SMILES strings to find similar molecules to.
        similarity: similarity threshold to use for the search. Value should be between 40 and 100. Defaults to 80.
        n_threads: Number of threads to use for searching the similar compounds. Defaults to 1.

    Returns:
        pd.DataFrame: a DataFrame with the similar molecules.
    """

    extracted = []

    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = []
        for smi in smiles:
            futures.append(executor.submit(get_similar_compounds, smi, similarity))
        for future in as_completed(futures):
            extracted.append(future.result())

    return pd.concat(extracted, ignore_index=True)


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
    assays_api = new_client.assay
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
            ).assign(variant_sequence=lambda x: x.variant_sequence.replace({None: "WT"}))
    else:
        activity_kwargs.pop("assay_chembl_id__in")
        raise ValueError(
            f"No assays found for the ids: {assay_chembl_ids} with the parameters: {activity_kwargs}"
        )
    return assays_df


def bioactivities_from_chembl(
    molecule_chembl_ids: Optional[list] = None,
    target_chembl_ids: Optional[list] = None,
    assay_chembl_ids: Optional[list] = None,
    document_chembl_ids: Optional[list] = None,
    **kwargs,
) -> pd.DataFrame:
    """Take a list of molecule chembl ids and get their respective bioactivities in ChEMBL.
    Args:
        molecule_chembl_id: list of molecule ChEMBL IDs to fecth bioactivities. Defaults to None.
        target_chembl_ids: list of target ChEMBL IDs to fetch bioactivities. Defaults to None.
        assay_chembl_ids: list of assay ChEMBL IDs to fetch bioactivities. Defaults to None.
        document_chembl_ids: list of document ChEMBL IDs to fetch bioactivities. Defaults to None.
        kwargs: example -> `standard_relation="=", assay_type__in=["B", "F"]`.
    Returns:
        pd.DataFrame: a DataFrame with the bioactivities.
    """
    activity_kwargs = {**kwargs}
    if molecule_chembl_ids is not None:
        activity_kwargs.update({"molecule_chembl_id__in": molecule_chembl_ids})
    if target_chembl_ids is not None:
        activity_kwargs.update({"target_chembl_id__in": target_chembl_ids})
    if assay_chembl_ids is not None:
        activity_kwargs.update({"assay_chembl_id__in": assay_chembl_ids})
    if document_chembl_ids is not None:
        activity_kwargs.update({"document_chembl_id__in": document_chembl_ids})
    activity_api = new_client.activity
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
