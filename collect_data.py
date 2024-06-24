import pandas as pd
import json
from pathlib import Path
from CompoundMapper.chembl import fetch_and_filter_workflow, molecule_info_from_chembl
from loguru import logger

all_dfs = []

with Path("jump_chemSmiSW_chembl_ids.json").open("r") as f:
    jump_SW_chembl_ids = json.load(f)

for key, values in jump_SW_chembl_ids.items():
    if values is None:
        continue
    print("Trying with ", key, values)
    try:
        chembl_data = fetch_and_filter_workflow(values, confidence_scores=[6, 7, 8, 9])
        jump_id_array = [key] * chembl_data.shape[0]
        unique_mol_chembl_ids = chembl_data["molecule_chembl_id"].unique().tolist()

        mol_data = molecule_info_from_chembl(unique_mol_chembl_ids)
        full_df = chembl_data.merge(mol_data, on="molecule_chembl_id", how="inner")
        full_df = full_df[
            ["molecule_chembl_id"]
            + [col for col in full_df if col != "molecule_chembl_id"]
        ]
        full_df.insert(0, "JUMP_ID", jump_id_array)

        if full_df.shape[0] != chembl_data.shape[0]:
            logger.warning(
                f"Merging went wrong??? check {key} - {unique_mol_chembl_ids}"
            )
        if (
            full_df["molecule_chembl_id"].nunique()
            != chembl_data["molecule_chembl_id"].nunique()
        ):
            logger.warning(
                f"Fetching the data went wrong??? check {key} - {unique_mol_chembl_ids}"
            )
        all_dfs.append(full_df)

    except ValueError:
        print(f"Error with {key}")
        continue

fulldf = pd.concat(all_dfs, ignore_index=True)
fulldf.to_csv("jump_chemSmiSW_chembl_data.csv", index=False)
