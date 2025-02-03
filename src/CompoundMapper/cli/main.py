"""Module containing the command line interface to get data from ChEMBL"""

import argparse
import json
from pathlib import Path

from .. import __version__
from ..logger import setup_logger
from .chembl_data_pipeline import aggregate_data, get_standardize_and_clean_workflow

DEFAULTS = {
    "molecule_ids": [],
    "target_ids": [],
    "assay_ids": [],
    "document_ids": [],
    "calculate_pchembl": False,  # store_true, default=False
    "output_path": "chembl_data.csv",
    "confidence_scores": [7, 8, 9],
    "bioactivity_type": ["Potency", "Kd", "Ki", "IC50", "AC50", "EC50"],
    "chirality": False,  # store_true, default=False
    "standard_relation": ["="],
    "assay_types": ["B", "F"],
    "log_level": "info",
    "chembl_version": None,
    "no_document_info": False,  # store_true, default=True
    "metadata_columns": [],
    "id_columns": [],
    "skip_not_aggregated": False,  # not(store_true), default=False
    "aggregate_mutants": False,  # not(store_true), default=False
    "save_recipe": True,
}

STORE_TRUE_ARGS = [
    "calculate_pchembl",
    "chirality",
    "no_document_info",
    "save_not_aggregated",
    "aggregate_mutants",
]


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
        default=DEFAULTS["output_path"],
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
        default=DEFAULTS["confidence_scores"],
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
        default=DEFAULTS["bioactivity_type"],
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
        default=DEFAULTS["standard_relation"],
        help=("Filter only bioactivities with the specified relation. Not implemented yet."),
        type=str,
    )
    parser.add_argument(
        "-at",
        "--assay_types",
        nargs="*",
        default=DEFAULTS["assay_types"],
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
        default=DEFAULTS["log_level"],
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
        default=DEFAULTS["chembl_version"],
    )
    parser.add_argument(
        "-no_doc",
        "--no_document_info",
        action="store_true",
        help=(
            "If passed, document information won't be included in the retrieved dataset. For example, "
            "year metadata will be missing, but requires one less API call. Defaults to False."
        ),
    )

    parser.add_argument(
        "-mcols",
        "--metadata_columns",
        nargs="*",
        default=DEFAULTS["metadata_columns"],
        help="Extra metadata columns to keep in the final dataframe, aggregated by ';'. Defaults to [].",
        type=str,
    )
    parser.add_argument(
        "-idcols",
        "--id_columns",
        nargs="*",
        default=DEFAULTS["id_columns"],
        help=(
            "Extra columns to use as identifiers for the aggregation. Passing `assay_chembl_id` to this "
            "argument, for example, will only aggregate the data if the compound is the same and the assay "
            "is the same. Saved data will still contain the original data separated by ';'. Defaults to []."
        ),
        type=str,
    )
    parser.add_argument(
        "-skip_not_agg",
        "--skip_not_aggregated",
        action="store_true",
        help="Skips saving the data before aggregation of same-molecule datapoint takes place.",
    )
    parser.add_argument(
        "-mutagg",
        "--aggregate_mutants",
        action="store_true",
        help=(
            "Aggregate data on targets regardless of the variant_sequence flag in ChEMBL, treating "
            "mutants as the same target. Regardless of the configuration, mutation data will be "
            "stored under `variant_sequence`. Default is False."
        ),
    )
    parser.add_argument(
        "-rec",
        "--save_recipe",
        help=(
            "Saves a json file with the parameters used to fetch the data. Useful for reproducibility. "
            "The file will be saved with the asme output path, but with the `_recipe.json` suffix."
            "Defaults to True."
        ),
        default=True,
        type=bool,
    )

    return parser.parse_args()


def main(args: argparse.Namespace) -> None:

    if args.standard_relation != ["="]:
        raise NotImplementedError("Fetching data using different relation types isn't implemented yet.")

    setup_logger(level=args.log_level.upper())
    output_path = Path(args.output_path)
    if not output_path.parent.exists():
        output_path.mkdir()

    df = get_standardize_and_clean_workflow(
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
        save_not_aggregated=(not args.skip_not_aggregated),
        add_document_info=(not args.no_document_info),
    )

    df = aggregate_data(
        df=df,
        chirality=args.chirality,
        chembl_version=args.chembl_version,
        metadata_cols=args.metadata_columns,
        extra_id_cols=args.id_columns,
        aggregate_mutants=args.aggregate_mutants,
        output_path=args.output_path,
    )

    # Save the recipe
    if args.save_recipe:
        output_name = output_path.stem
        recipe_path = output_path.parent / f"{output_name}_recipe.json"

        configs = vars(args)
        command_vals = []
        for k, v in configs.items():
            if isinstance(v, bool):
                if DEFAULTS[k] != v:
                    command_vals.append((f"--{k}" if k in STORE_TRUE_ARGS else f"--{k} {v}"))
            elif isinstance(v, list):
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{k} {' '.join([str(i) for i in v])}")
            elif isinstance(v, str):
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{k} {v}")

        command = "getchembl " + " ".join(command_vals)
        configs.update({"CompoundMapper version": __version__})
        configs.update({"command": command})

        with open(recipe_path, "w") as f:
            json.dump(configs, f, indent=2)

    return df


def main_exe() -> None:
    args = parse_arguments()
    main(args)
