"""Module containing the command line interface to get data from ChEMBL"""

import argparse
from pathlib import Path

import numpy as np
from chemFilters.chem.standardizers import ChemStandardizer
from UniProtMapper import ProtMapper

from ..chembl.processing import fetch_and_filter_workflow
from ..core.smiles_utils import clean_mixtures
from ..core.stats_make import process_repeat_mols, repeated_indices_from_array_series
from ..logger import logger, setup_logger
from .fp_utils import calculate_mixed_FPs
from .workflow import aggregate_data, fetch_standardize_and_clean_workflow


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
        help="ChEMBL molecule IDs to download data from.",
        type=str,
    )
    parser.add_argument(
        "-tids",
        "--target_ids",
        dest="target_ids",
        nargs="*",
        required=False,
        help="ChEMBL target IDs to retrieve data from.",
        type=str,
    )
    parser.add_argument(
        "-asids",
        "--assay_ids",
        dest="assay_ids",
        nargs="*",
        required=False,
        help="ChEMBL assay IDs to retrieve data from.",
        type=str,
    )
    parser.add_argument(
        "-dids",
        "--document_ids",
        dest="document_ids",
        nargs="*",
        required=False,
        help="ChEMBL document IDs to download data from.",
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
    parser.add_argument(
        "-v",
        "--chembl_version",
        type=int,
        help="chembl_version: specify latest ChEMBL release to extract data from (e.g., 28). Defaults to None.",
        default=None,
    )
    parser.add_argument(
        "-mcols",
        "--metadata_cols",
        nargs="*",
        default=[],
        help="Extra metadata columns to keep in the final dataframe, aggregated by ';'. Defaults to [].",
        type=str,
    )
    parser.add_argument(
        "-noaggr",
        "--save_not_aggregated",
        action="store_true",
        help="Save the data before aggregating the repeated molecules.",
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:

    if args.standard_relation != ["="]:
        raise NotImplementedError("Fetching data using different relation types isn't implemented yet.")

    setup_logger(level=args.log_level.upper())
    output_path = Path(args.output_path)
    if not output_path.parent.exists():
        output_path.mkdir()

    df = fetch_standardize_and_clean_workflow(
        molecule_ids=args.molecule_ids,
        target_ids=args.target_ids,
        assay_ids=args.assay_ids,
        document_ids=args.document_ids,
        calculate_pchembl=args.calculate_pchembl,
        output_path=args.output_path,
        confidence_scores=args.confidence_scores,
        bioactivity_type=args.bioactivity_type,
        standard_relation=args.standard_relation,
        assay_types=args.assay_types,
        chembl_version=args.chembl_version,
        save_not_aggregated=args.save_not_aggregated,
    )

    df = aggregate_data(
        df=df,
        chirality=args.chirality,
        chembl_version=args.chembl_version,
        metadata_cols=args.metadata_cols,
        output_path=args.output_path,
    )
    return df


def main_exe() -> None:
    args = parse_arguments()
    main(args)
