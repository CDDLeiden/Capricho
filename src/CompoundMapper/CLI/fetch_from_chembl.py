"""Module containing the command line interface to get data from ChEMBL"""

import argparse
from pathlib import Path

import numpy as np
from chemFilters.chem.standardizers import ChemStandardizer
from UniProtMapper import ProtMapper

from ..chembl import fetch_and_filter_workflow
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import process_repeat_mols, repeated_indices_from_array_series
from ..logger import logger, setup_logger
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
        dest="confidence_scores",
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
    parser.add_argument(
        "-chiral",
        "--chirality",
        dest="chirality",
        help=(
            "Consider chirality when calculating the fingerprints. Useful if you want to "
            "differentiate between enantiomers while aggregating the data. Default is False."
        ),
        action="store_true",
    )
    parser.add_argument(
        "-rel",
        "--standard_relation",
        nargs="*",
        default=["="],
        help=("Filter only bioactivities with the specified relation. Not implemented yet."),
        type=str,
    )
    parser.add_argument(
        "-at",
        "--assay_types",
        nargs="*",
        default=["B", "F"],
        help=(
            "Assay types to filter. Defaults to both functional (F) and binding assays (B). "
            "Other assay types: ADME (A), Toxicity (T) and Physichochemical (P)."
        ),
        type=str,
    )
    parser.add_argument(
        "-log",
        "--log-level",
        dest="log_level",
        default="info",
        choices=["info", "debug", "warning", "error", "critical"],
        help=(
            "Set the logging level. Defaults to info. "
            "Choose between: info, debug, warning, error, critical."
        ),
        type=str,
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:

    if args.standard_relation != ["="]:
        raise NotImplementedError("Fetching data using different relation types isn't implemented yet.")

    setup_logger(level=args.log_level.upper())
    output_path = Path(args.output_path)
    if not output_path.parent.exists():
        output_path.mkdir()

    # since we work with pchembl values, standard values reported as -pXC50, -Log XC50, etc. will be renamed
    bioactivity_type_rename_dict = {
        **{f"p{bio}": bio for bio in args.bioactivity_type},
        **{f"Log {bio}": bio for bio in args.bioactivity_type},
        **{f"-Log {bio}": bio for bio in args.bioactivity_type},
    }

    full_df = fetch_and_filter_workflow(
        molecule_chembl_ids=args.molecule_ids,
        target_chembl_ids=args.target_ids,
        confidence_scores=args.confidence_scores,
        assay_types=args.assay_types,
        calculate_pchembl=args.calculate_pchembl,
    )
    # full_df.to_csv(output_path.with_suffix(".begin.csv"), index=False)

    stdzer = ChemStandardizer(from_smi=True, n_jobs=8, verbose=False)
    queried_df = (
        full_df.assign(  # rename bioactivities & filter by preferred bioactivity type
            standard_type=lambda x: x["standard_type"].replace(bioactivity_type_rename_dict)
        )
        .query("standard_type.isin(@args.bioactivity_type)")
        # standardize the smiles & clean possible solvents & salts from the string
        .assign(standard_smiles=lambda x: stdzer(x["canonical_smiles"]))
        .dropna(subset=["standard_smiles"])  # drop if no structure is found
        .query("standard_smiles.notna()")
        .assign(final_smiles=lambda x: x["standard_smiles"].apply(clean_mixtures))
        .drop(columns="standard_smiles")
        .rename(columns={"final_smiles": "standard_smiles"})
        # Drop cols that won't be used; we'll use `standard_<colname>` instead
        .drop(
            columns=[
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
        .reset_index(drop=True)
        .copy()
    )
    # drop rows with missing pchembl values
    no_pchembl_idxs = queried_df.query("pchembl_value.isna()").index
    logger.info(f"Dropping {len(no_pchembl_idxs)} rows missing pchembl values.")
    queried_df = queried_df.drop(index=no_pchembl_idxs).reset_index(drop=True)

    # find remaining mixtures in the data
    mask = queried_df["standard_smiles"].str.contains(".", regex=False)
    n_mixtures = mask.sum()
    if n_mixtures > 0:
        logger.info(f"Number of mixtures: {mask.sum()}")
        mixture_idxs = np.where(mask)[0]  # drop data points where smiles are mixtures
        queried_df = queried_df.drop(index=mixture_idxs).reset_index(drop=True)

    # calculate the fingerprints to identify repeat molecules & aggregate data accordingly
    fps = calculate_mixed_FPs(
        queried_df["standard_smiles"].tolist(), n_jobs=4, morgan_kwargs={"useChirality": args.chirality}
    )
    # queried_df.to_csv(output_path.with_suffix(".midway.csv"), index=False)
    queried_df = queried_df.assign(fps=fps)
    repeats_idxs = repeated_indices_from_array_series(queried_df["fps"])
    final_data = process_repeat_mols(queried_df, repeats_idxs, solve_strat="keep", chirality=args.chirality)
    final_data.to_csv(output_path, index=False)
    return final_data  # could be used as a function


def main_exe() -> None:
    args = parse_arguments()
    main(args)
