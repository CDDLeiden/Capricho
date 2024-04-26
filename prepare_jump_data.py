from chemFilters.chem.standardizers import ChemStandardizer

import pandas as pd
from loguru import logger
from rdkit import Chem

stdzer_iso = ChemStandardizer(method='chembl', isomeric=True, from_smi=True, n_jobs=8)
stdzer_not_iso = ChemStandardizer(method='chembl', isomeric=False, from_smi=True, n_jobs=8)

fname = "JUMPCP_compounds_processed.csv.gz"
meta_col = "Metadata_JCP2022"
inchikey_col = "Standardized_inchikey"


def inchikey_from_smi(smi):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smi))


def unichem_res_to_df(res: dict):
    df = pd.DataFrame.from_dict(res)
    comparison_cols = [c for c in df.columns if c.startswith("comparison_")]
    df = df.assign(
        exact_match=df[comparison_cols].apply(lambda x: x.all(), axis=1)
    ).query(
        "id == 1"  # corresponds to ChEMBL
    )
    return df


if inchikey_col not in pd.read_csv(fname, nrows=1).columns:
    subset = pd.read_csv(fname).assign(
        Standardized_inchikey=lambda x: x["Standardized_SMILES"].apply(
            inchikey_from_smi
        )
    )
    subset.to_csv(fname, index=False)
else:
    subset = pd.read_csv(fname)
    
if all(['chembl_iso_smiles' not in subset.columns, 'chembl_not_iso_smiles' not in subset.columns]):
    subset = subset.assign(
        chembl_iso_smiles=lambda x: stdzer_iso(x["smiles"]),
        chembl_not_iso_smiles=lambda x: stdzer_not_iso(x["smiles"])
    )
    subset.to_csv(fname, index=False)

# assert all inchikeys and JUMP identifiers are unique
if subset[inchikey_col].nunique() != subset.shape[0]:
    logger.warning("There are possibly duplicated compounds!")
if subset[meta_col].nunique() != subset.shape[0]:
    logger.warning("There are possibly duplicated IDs!")

mol_cols = [
    "hierarchy_active_id",
    "hierarchy_molecule_id",
    "hierarchy_parent_id",
]

subset.head(1)