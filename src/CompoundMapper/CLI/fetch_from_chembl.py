"""Module containing the command line interface to get data from ChEMBL"""

import argparse
from pathlib import Path

import numpy as np
from chemFilters.chem.standardizers import ChemStandardizer
from UniProtMapper import ProtMapper

from CompoundMapper.core.pandas_helper import assign_stats

from ..chembl import fetch_and_filter_workflow, molecule_info_from_chembl
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import process_repeat_mols, repeated_indices_from_array_series
from ..logger import logger
from .fp_utils import calculate_mixed_FPs


def fetch_names(chembl_ids: str):
    retriever = ProtMapper()
    results, failed = retriever(chembl_ids, from_db="ChEMBL", fields=["protein_name", "organism_name"])
    return results, failed


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="ChEMBL fetcher",
        description="A command line interface to filter, download and process data from ChEMBL.",
    )
    parser.add_argument(
        "-mids",
        "--molecule_ids",
        dest="molecule_ids",
        nargs="*",
        required=False,
        help="Chembl molecule IDs to download data from.",
        type=str,
    )
    parser.add_argument(
        "-tids",
        "--target_ids",
        dest="target_ids",
        nargs="*",
        required=False,
        help="UniProt target IDs to download data from.",
        type=str,
    )
    parser.add_argument(
        "-calc",
        "--calculate_pchembl",
        action="store_true",
        help=(
            "Calculate the pChEMBL (pXC50) values when not reported for bioactivities "
            "reported in nM, µM or uM. Default is False."
        ),
    )
    parser.add_argument(
        "-o",
        "--output_path",
        type=str,
        help=(
            "Path to save the output files. If not provided, "
            "will create a folder named 'chembl_data' in the current directory."
        ),
        default="chembl_data",
    )
    parser.add_argument(
        "-c",
        "--confidence_scores",
        nargs="*",
        required=False,
        help=(
            "Confidence scores to filter the bioactivities. "
            "If not provided, will fetch only score 9 bioactivities."
        ),
        default=[8, 9],
        type=int,
    )
    parser.add_argument(
        "-biotype",
        "--bioactivity_type",
        dest="bioactivity_type",
        help=(
            "Type of bioactivity to filter. If left empty, will fetch for all types. "
            "Examples of bioactivity types: `Potency`, `Kd`, `Ki`, `IC50`, `AC50`, `EC50`."
        ),
        default=["Potency", "Kd", "Ki", "IC50", "AC50", "EC50"],
        nargs="*",
        type=str,
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:

    output_path = Path(args.output_path)
    if not output_path.parent.exists():
        output_path.mkdir()

    # since we work with pchembl values, standard values reported as -pXC50, -Log XC50, etc. will be renamed
    bioactivity_type_rename_dict = {
        **{f"p{bio}" for bio in args.bioactivity_type},
        **{f"Log {bio}": bio for bio in args.bioactivity_type},
        **{f"-Log {bio}": bio for bio in args.bioactivity_type},
    }

    chembl_data = fetch_and_filter_workflow(
        molecule_chembl_ids=args.molecule_ids,
        target_chembl_ids=args.target_ids,
        confidence_scores=args.confidences,
    )

    unique_mol_chembl_ids = chembl_data["molecule_chembl_id"].unique().tolist()
    mol_data = molecule_info_from_chembl(unique_mol_chembl_ids)
    full_df = chembl_data.merge(mol_data, on="molecule_chembl_id", how="inner")
    # reorder the columns...
    full_df = full_df[["molecule_chembl_id"] + [col for col in full_df if col != "molecule_chembl_id"]]

    stdzer = ChemStandardizer(from_smi=True, n_jobs=8, verbose=False)

    queried_df = (
        full_df.copy()  # rename bioactivities & filter by preferred bioactivity type
        .assign(standard_type=lambda x: x["standard_type"].replace(bioactivity_type_rename_dict))
        .query("standard_type.isin(@args.bioactivity_type)")
        .assign(  # standardize the smiles & clean possible solvents & salts from the string
            standard_smiles=lambda x: stdzer(x["canonical_smiles"]),
            final_smiles=lambda x: x["standard_smiles"].apply(clean_mixtures),
        )
        .drop(columns="standard_smiles")
        .rename(columns={"final_smiles": "standard_smiles"})
        .drop(
            columns=[  # columns that won't be used; we'll use `standard_<colname>` instead
                "type",
                "relation",
                "units",
                "value",
                "standard_value",  # we'll use pchembl instead
                "type",
                "canonical_smiles",
                # "description",
            ]
        )
    )
    # drop rows with missing pchembl values
    no_pchembl_idxs = queried_df.query("pchembl_value.isna()").index
    logger.info(f"Dropping {len(no_pchembl_idxs)} rows missing pchembl values.")
    queried_df = queried_df.drop(index=no_pchembl_idxs).reset_index(drop=True)

    mask = queried_df["standard_smiles"].str.contains(".", regex=False)
    logger.info(f"Number of mixtures: {mask.sum()}")
    mixture_idxs = np.where(mask)[0]  # drop data points where smiles are mixtures
    queried_df = queried_df.drop(index=mixture_idxs).reset_index(drop=True)

    # calculate the fingerprints to identify repeat molecules & aggregate data accordingly
    mol_fp = calculate_mixed_FPs(queried_df["standard_smiles"].tolist(), n_jobs=4)
    repeats_idxs = repeated_indices_from_array_series(mol_fp)
    final_data, final_cols = process_repeat_mols(queried_df, repeats_idxs, solve_strat="keep")
    final_data = final_data[final_cols].drop_duplicates().reset_index(drop=True)
    final_data = assign_stats(final_data)

    final_data.to_csv(output_path, index=False)


def main_exe() -> None:
    args = parse_arguments()
    main(args)
