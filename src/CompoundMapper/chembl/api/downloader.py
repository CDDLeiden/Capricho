"""Module holding functionalities for the ChEMBL API using chembl downloader as the backend."""

import json
from pathlib import Path
from textwrap import dedent
from typing import List, Optional, Sequence, Tuple, Union

import pandas as pd
import pystow

from chembl_downloader import download_extract_sqlite, latest, query
from chembl_downloader.api import _find_sqlite_file

from ...logger import logger

PYSTOW_PARTS = ["chembl"]
PYSTOW_CONFIG = {"name": "chembl_downloader_config_{version}.json"}


def _get_config_file(version: Optional[Union[int, str]] = None) -> Path:
    version = version if version is not None else latest()
    return pystow.join(*(PYSTOW_PARTS), PYSTOW_CONFIG["name"].format(version=version))


def check_and_download_chembl_db(
    prefix: Optional[Sequence[str]] = None,
    version: Optional[Union[int, str]] = None,
) -> pd.DataFrame:
    """Check if the ChEMBL database is present and download it if not.

    Args:
        prefix: Optional prefix for an alternative data directory. If passed, will create
            a new configuration file under `~/.data/chembl_downloader_config_{version}.json`
            pointing to the new data directory. Defaults to None.
        version: Optional ChEMBL version to download. If not provided, will download the latest
            available version. Defaults to None.

    Returns:
        Path to the ChEMBL SQLite database
    """
    # if present, config file override the default path, unless a prefix is defined
    version = version if version is not None else latest()
    config_file = _get_config_file(version)
    configs = {"prefix": (prefix if prefix is not None else prefix), "version": version}

    if config_file.exists():  # only exists if that version was downloaded to custom path before
        logger.info(f"Loading ChEMBL configuration from:\n\t{config_file}")
        configs = json.loads(config_file.read_text())
        logger.debug(f"Loaded configuration:\n{json.dumps(configs, indent=2)}")

    sql_path = _find_sqlite_file(pystow.join(*(configs["prefix"]), f"{configs['version']}/data"))
    if not sql_path.exists():
        logger.info(f"Downloading and extracting ChEMBL version {version} into:\n\t{sql_path}")
        download_extract_sqlite(**configs)
        if prefix is not None:
            config_file.write_text(json.dumps(configs, indent=2))
    else:
        logger.info(f"Loading local ChEMBL database at:\n\t{sql_path}")

    return configs


def get_full_activity_data_sql(
    molecule_chembl_ids: Optional[Union[list, str]] = None,
    target_chembl_ids: Optional[Union[list, str]] = None,
    assay_chembl_ids: Optional[Union[list, str]] = None,
    document_chembl_ids: Optional[Union[list, str]] = None,
    standard_relation: Optional[List[str]] = None,
    standard_type: Optional[List[str]] = None,
    confidence_scores: Union[list, Tuple] = (9, 8),
    assay_types: Union[list, Tuple] = ("B", "F"),
    chembl_version: Optional[int] = None,
    additional_fields: Optional[List[str]] = None,
    prefix: Optional[Sequence[str]] = None,
    version: Optional[Union[int, str]] = None,
) -> pd.DataFrame:
    """Retrieve ChEMBL bioactivity data from any combination of molecule, target, assay, or document IDs.
    Data is retrieved using the ChEMBL downloader. Merges are performed on the SQL query level and a
    DataFrame is returned with the bioactivity data.

    Args:
        molecule_chembl_ids: list of ChEMBL molecule IDs to fetch data for. Defaults to None.
        target_chembl_ids: list of ChEMBL target IDs to fetch data for. Defaults to None.
        assay_chembl_ids: list of ChEMBL assay IDs to fetch data for. Defaults to None.
        document_chembl_ids: list of ChEMBL document IDs to fetch data for. Defaults to None.
        standard_relation: Optional filter for standard relation types (e.g., ["=", "<", ">"])
        standard_type: Optional filter for activity types (e.g., ["IC50", "Ki", "EC50"])
        confidence_scores: list of confidence scores to filter the fetched assay data.
            Defaults to (9, 8).
        assay_types: list of assay types to be fetched from ChEMBL. Defaults to binding (B) and
            functional (F) data.
        chembl_version: Not to confuse for `version`. This is the ChEMBL release number used to
            filter the data. Defaults to None.
        additional_fields: Optional list of additional fields to include in the sql query. E.g.:
            ["vs.sequence"], to retrieve the sequence of the variant, if available. Defaults to None.
        prefix: Optional prefix for an alternative data directory. If passed, will create
            a new configuration file under `~/.data/chembl_downloader_config_{version}.json`
            pointing to the new data directory. Defaults to None.
        version: ChEMBL database to be downloaded and used by ChEMBL downloader. If not provided,
            will download the latest available version. Defaults to None.

    Returns:
        pd.DataFrame: a DataFrame with the bioactivity data.
    """
    # Download the ChEMBL database if not present
    downloader_configs = check_and_download_chembl_db(prefix=prefix, version=version)
    where_conditions = []  # store the WHERE conditions

    # handle optional filters for molecule, target, assay, and document IDs
    if molecule_chembl_ids is not None:
        if isinstance(molecule_chembl_ids, str):
            molecule_ids = [molecule_chembl_ids]
        else:
            molecule_ids = molecule_chembl_ids
        m_placeholders = ", ".join([f"'{id}'" for id in molecule_ids])
        where_conditions.append(f"md.chembl_id IN ({m_placeholders})")

    if target_chembl_ids is not None:
        if isinstance(target_chembl_ids, str):
            target_ids = [target_chembl_ids]
        else:
            target_ids = target_chembl_ids
        t_placeholders = ", ".join([f"'{id}'" for id in target_ids])
        where_conditions.append(f"td.chembl_id IN ({t_placeholders})")

    if document_chembl_ids is not None:
        if isinstance(document_chembl_ids, str):
            document_ids = [document_chembl_ids]
        else:
            document_ids = document_chembl_ids
        d_placeholders = ", ".join([f"'{id}'" for id in document_ids])
        where_conditions.append(f"d.chembl_id IN ({d_placeholders})")

    if assay_chembl_ids is not None:
        if isinstance(assay_chembl_ids, str):
            assay_ids = [assay_chembl_ids]
        else:
            assay_ids = assay_chembl_ids
        a_placeholders = ", ".join([f"'{id}'" for id in assay_ids])
        where_conditions.append(f"a.chembl_id IN ({a_placeholders})")

    # Build field list
    base_fields = [
        "act.activity_id AS activity_chembl_id",
        "a.chembl_id AS assay_chembl_id",
        "a.description AS assay_description",
        "a.assay_type",
        "a.assay_organism",
        "a.assay_category",
        "a.assay_tax_id",
        "a.assay_strain",
        "a.assay_tissue",
        "a.assay_cell_type",
        "a.assay_subcellular_fraction",
        "a.bao_format",
        "a.confidence_score",
        "md.chembl_id AS molecule_chembl_id",
        "md.first_in_class",
        "md.chirality",
        "md.oral",
        "md.prodrug",
        "md.max_phase",
        "md.therapeutic_flag",
        "md.withdrawn_flag",
        "md.indication_class",
        "act.standard_flag",
        "act.standard_relation",
        "act.standard_type",
        "act.standard_units",
        "act.standard_value",
        "act.pchembl_value",
        "td.chembl_id AS target_chembl_id",
        "td.organism",
        "cs.canonical_smiles",
        "cs.standard_inchi_key",
        # "cs.molregno",
        # "cs.active_molregno",
        # "cs.parent_molregno",
        "act.data_validity_comment",
        "act.potential_duplicate",
        "act.data_validity_comment",
        "d.chembl_id AS document_chembl_id",
        "d.doc_type",
        "d.authors",
        "d.doi",
        "d.journal",
        "d.volume",
        "d.year",
        "d.title",
        "vs.mutation",
        "d.chembl_release_id",
    ]

    if additional_fields:
        all_fields = base_fields + additional_fields
    else:
        all_fields = base_fields

    fields_clause = ",\n            ".join(all_fields)

    where_conditions.append("act.standard_value IS NOT NULL")
    where_conditions.append("act.standard_relation IS NOT NULL")

    if standard_relation:
        relation_placeholders = ", ".join([f"'{rel}'" for rel in standard_relation])
        where_conditions.append(f"act.standard_relation IN ({relation_placeholders})")

    if standard_type:
        type_placeholders = ", ".join([f"'{type}'" for type in standard_type])
        where_conditions.append(f"act.standard_type IN ({type_placeholders})")

    if confidence_scores:
        confidence_placeholders = ", ".join([f"{score}" for score in confidence_scores])
        where_conditions.append(f"a.confidence_score IN ({confidence_placeholders})")

    if assay_types:
        atype_placeholders = ", ".join([f"'{type}'" for type in assay_types])
        where_conditions.append(f"a.assay_type IN ({atype_placeholders})")

    if chembl_version:
        where_conditions.append(f"(d.chembl_release_id IS NULL OR d.chembl_release_id <= {chembl_version})")

    where_clause = " AND\n            ".join(where_conditions)

    query_str = dedent(  # Build complete query
        f"""\
        SELECT
            {fields_clause}
        FROM molecule_dictionary md
        JOIN compound_structures cs ON md.molregno = cs.molregno
        JOIN activities act ON md.molregno = act.molregno
        JOIN docs d ON act.doc_id = d.doc_id
        JOIN assays a ON act.assay_id = a.assay_id
        LEFT JOIN variant_sequences vs ON a.variant_id = vs.variant_id
        JOIN target_dictionary td ON a.tid = td.tid
        WHERE
            {where_clause}
        ORDER BY
            md.chembl_id, act.standard_type
    """
    )

    logger.debug(f"Generated SQL query:\n{query_str}")

    return query(
        query_str,
        version=downloader_configs["version"],
        prefix=downloader_configs["prefix"],
    )
