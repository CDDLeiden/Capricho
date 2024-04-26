import json
import time
import pandas as pd
from loguru import logger
from rdkit import Chem
from tqdm import tqdm
from smallworld_api import SmallWorld

fname = "JUMPCP_compounds_processed.csv.gz"
meta_col = "Metadata_JCP2022"
smiles_col = "chembl_not_iso_smiles"

print('Using parameters:')
print({'fname': fname, 'meta_col': meta_col, 'smiles_col': smiles_col})


def inchikey_from_smi(smi):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smi))


if smiles_col not in pd.read_csv(fname, nrows=1).columns:
    subset = pd.read_csv(fname).assign(
        Standardized_inchikey=lambda x: x["Standardized_SMILES"].apply(
            inchikey_from_smi
        )
    )
    subset.to_csv(fname, index=False)
else:
    subset = pd.read_csv(fname)

# assert all inchikeys and JUMP identifiers are unique
if subset[smiles_col].nunique() != subset.shape[0]:
    logger.warning("There are possibly duplicated compounds!")
if subset[meta_col].nunique() != subset.shape[0]:
    logger.warning("There are possibly duplicated IDs!")

mol_cols = [
    "hierarchy_active_id",
    "hierarchy_molecule_id",
    "hierarchy_parent_id",
]

sw = SmallWorld()
db = sw.chembl_dataset
sw.configure_logger(level="WARNING")

default_submission = {
    "sdist": 12,  # maximum scored distance; will limit the search according to the `scores`
    "dist": 8,  # maximum anonymous (topology) distance; will limit the search to [tdn, tup, rdn, rup, ldn, lup]=dist
    "tdn": 6,  # number of terminal down (tdn) edits
    "tup": 6,  # number of terminal up (tup) edits
    "rdn": 6,  # number of ring down (rdn) edits
    "rup": 2,  # number of ring up (rup) edits
    "ldn": 2,  # number of linker down (ldn) edits
    "lup": 2,  # number of linker up (lup) edits
    "maj": 6,  # number of major element transmutations
    "min": 6,  # number of minor element transmutations
    "sub": 6,  # number of substitutions (heavy degree change)
    "scores": "Atom Alignment,ECFP4,Daylight",
}

search_params = {
    "db": db,
    "dist": 0,  # This limits the search to sum([tdn, tup, rdn, rup, ldn, lup]) = 0
    "sdist": 1,  # defined restrictions already stric enough so we can allow scoring distances to be greater
    "length": 30,  # results will be always less than 30
    "maj": 0, #  no major element transmutations
    "min": 0, #  no minor element transmutations
    "sub": 0, #  no heavy degree changes (Hydrogen substitution)
}

search_params = {**sw.default_submission, **search_params}

print('SW parameters:')
print(search_params)

results_dict = {}
smiles_dict = {}
time_seconds = []
for jumpID, smi in tqdm(
    zip(subset[meta_col], subset[smiles_col]), total=subset.shape[0]
):
    start = time.time()
    results: pd.DataFrame = sw.search(smi, **search_params)
    end = time.time()
    time_seconds.append(end - start)
    
    if results.empty:
        results_dict.update({jumpID: None})
        smiles_dict.update({jumpID: None})
    else:
        results_dict.update({jumpID: results['id'].tolist()})
        smiles_dict.update({jumpID: results['hitSmiles'].tolist()})
    time.sleep(0.25)

with open("jump_chemSmiSW_chembl_ids.json", "w") as json_file:
    json.dump(results_dict, json_file, indent=4)
    
with open("jump_chemSmiSW_smiles.json", "w") as json_file:
    json.dump(smiles_dict, json_file, indent=4)
    
# save the time it took for each request as pandas series
time_series = pd.Series(time_seconds, name="time_seconds")
time_series.to_csv("jump_SW_time.csv", index=False)
