"""Parsing methods to extract information retrieved from the ChEMBL API."""

from ..logger import logger


def parse_compound_response(r: dict, compound_id: str) -> dict:
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
