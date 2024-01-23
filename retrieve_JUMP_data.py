import json
import time

import numpy as np
import pandas as pd
from loguru import logger
from rdkit import Chem
from tqdm import tqdm

from CompoundMapper.chembl import molecule_info_from_chembl
from CompoundMapper.unichem import UniChem

unichem = UniChem()
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

results_dict = {}
counter = 0
for jumpID, inchikey in tqdm(
    zip(subset[meta_col], subset[inchikey_col]), total=subset.shape[0]
):
    counter += 1
    if counter == 10:
        time.sleep(1.5)
        counter = 0

    res = unichem.get_connectivity(inchikey)
    if res is not None:
        uchem_df = unichem_res_to_df(res)
        if uchem_df.empty:
            # logger.info(f"No ChEMBL molecules found for {jumpID}")
            results_dict.update({jumpID: None})
            continue
    else:
        # logger.info(f"No molecules found for {jumpID}")
        results_dict.update({jumpID: None})
        continue
    unichem_chembl_ids = uchem_df["compoundId"].tolist()
    try:
        mol_data = molecule_info_from_chembl(unichem_chembl_ids)
    # if it finds a None:
    except AttributeError:
        logger.warning(f"ERROR FOR ID {jumpID}")
        results_dict.update({jumpID: None})
        continue
    all_mol_ids = np.unique(
        mol_data[mol_cols].values.ravel().tolist() + unichem_chembl_ids
    ).tolist()
    results_dict.update({jumpID: all_mol_ids})

# save results_dict as a json
with open("jump_chembl_ids.json", "w") as json_file:
    json.dump(results_dict, json_file, indent=4)
