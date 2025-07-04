"""Module containing the command line interface to get data from ChEMBL"""

import argparse
import json
import sys
from pathlib import Path

from chembl_downloader import latest

from .. import __version__
from ..chembl.api.downloader import check_and_download_chembl_db
from ..chembl.api.sql_explorer import explorer_main
from ..core.default_fields import (
    DATA_DROPPING_COMMENT,
    DATA_PROCESSING_COMMENT,
    DEFAULT_ASSAY_MATCH_FIELDS,
)
from ..logger import logger, setup_logger
from .chembl_data_pipeline import aggregate_data, get_standardize_and_clean_workflow

DEFAULTS = {
    "molecule_ids": [],
    "target_ids": [],
    "assay_ids": [],
    "document_ids": [],
    "calculate_pchembl": False,
    "output_path": "chembl_data.csv",
    "confidence_scores": [7, 8, 9],
    "bioactivity_type": ["Potency", "Kd", "Ki", "IC50", "AC50", "EC50"],
    "chirality": False,
    "standard_relation": ["="],
    "assay_types": ["B", "F"],
    "log_level": "info",
    "chembl_release": None,
    "metadata_columns": [],
    "id_columns": [],
    "skip_not_aggregated": False,
    "aggregate_mutants": False,
    "skip_recipe": False,
    "drop_unassigned_chiral": False,
    "curate_annotation_errors": False,
    "chembl_backend": "downloader",
    "chembl_version": None,
    "save_dropped": False,
    "require_doc_date": False,
    "max_assay_size": None,
    "min_assay_size": None,
    "max_assay_match": False,
    "min_assay_overlap": 0,
    "strict_mutant_removal": False,
    "keep_flagged_data": False,
    "compound_equality": "connectivity",
}

STORE_TRUE_ARGS = [
    "calculate_pchembl",
    "chirality",
    "skip_not_aggregated",
    "aggregate_mutants",
    "skip_recipe",
    "drop_unassigned_chiral",
    "curate_annotation_errors",
    "save_dropped",
    "require_doc_date",
    "max_assay_match",
    "strict_mutant_removal",
    "keep_flagged_data",
]

STORE_FALSE_ARGS = ["skip_recipe"]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="CompoundMapper",
        description="A command line interface to filter, download and process data from ChEMBL.",
    )
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    # Explore the downloaded ChEMBL SQL database
    explore_parser = subparsers.add_parser(
        "explore",
        help=(
            "Explore the ChEMBL SQL database downloaded using chembl_downloader. For a visual "
            "inspection of the latest ChEMBL schema, consider checking: https://www.ebi.ac.uk/chembl/db_schema"
        ),
    )
    explore_parser.add_argument(
        "--version",
        "-v",
        required=False,
        help="ChEMBL version to use. If not provided, will use the latest version.",
        type=int,
    )
    explore_parser.add_argument(
        "--list-tables",
        "-list",
        dest="list_tables",
        action="store_true",
        help="List all tables within the SQL database and exit.",
    )
    explore_parser.add_argument(
        "--table",
        "-t",
        help="Explore a specific table",
    )
    explore_parser.add_argument(
        "--search-column",
        "-search",
        dest="search_column",
        help="Search for tables containing column name pattern.",
    )
    explore_parser.add_argument(
        "--query",
        "-q",
        help="Run a custom SQL query",
    )
    explore_parser.add_argument(
        "-log",
        "--log-level",
        dest="log_level",
        default=DEFAULTS["log_level"],
        choices=["trace", "debug", "info", "warning", "error", "critical"],
        help="Set the logging level. Defaults to info",
    )

    # Download ChEMBL SQL database using chembl_downloader
    download_parser = subparsers.add_parser(
        "download", help="Download ChEMBL SQL database using chembl_downloader."
    )
    download_parser.add_argument(
        "--version",
        "-v",
        dest="version",
        help="ChEMBL version download and extract. If not provided, will use the latest version.",
        default=None,
        type=int,
    )
    download_parser.add_argument(
        "--prefix",
        "-p",
        dest="prefix",
        help=(
            "Path to be used by pystow to store the data. If a custom prefix is passed, a config.json "
            "will created under `~/.data/` pointing to the custom path whenever the same ChEMBL version "
            "(--version) is referred to. Defaults to None, saving under `~/.data/chembl/."
        ),
        nargs="*",
        default=None,
        type=str,
    )
    download_parser.add_argument(
        "-log",
        "--log-level",
        dest="log_level",
        default=DEFAULTS["log_level"],
        choices=["trace", "debug", "info", "warning", "error", "critical"],
        help="Set the logging level. Defaults to info",
    )

    parser.add_argument(
        "-mids",
        "--molecule-ids",
        dest="molecule_ids",
        nargs="*",
        required=False,
        help="ChEMBL molecule IDs to download data from.",
        type=str,
    )
    parser.add_argument(
        "-tids",
        "--target-ids",
        dest="target_ids",
        nargs="*",
        required=False,
        help="ChEMBL target IDs to retrieve data from.",
        type=str,
    )
    parser.add_argument(
        "-asids",
        "--assay-ids",
        dest="assay_ids",
        nargs="*",
        required=False,
        help="ChEMBL assay IDs to retrieve data from.",
        type=str,
    )
    parser.add_argument(
        "-dids",
        "--document-ids",
        dest="document_ids",
        nargs="*",
        required=False,
        help="ChEMBL document IDs to download data from.",
        type=str,
    )
    parser.add_argument(
        "-calc",
        "--calculate-pchembl",
        action="store_true",
        help=(
            "Calculate the pChEMBL (pXC50) values when not reported for bioactivities "
            "reported in nM, µM or uM. Default is False."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-path",
        type=str,
        help=(
            "Path to save the output files. If not provided, "
            "will create a folder named 'chembl_data' in the current directory."
        ),
        default=DEFAULTS["output_path"],
    )
    parser.add_argument(
        "-c",
        "--confidence-scores",
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
        "--bioactivity-type",
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
        "-duchi",
        "--drop-unassigned-chiral",
        dest="drop_unassigned_chiral",
        help=(
            "Define the behavior for dealing with entries that have unassigned chiral centers. When passed, "
            "entries with unassigned chiral center will be dropped as part of the dataset curation. We recommend "
            "passing this flag if you're also passing the `chirality` flag. Default is False."
        ),
        action="store_true",
    )
    parser.add_argument(
        "-cae",
        "--curate-annotation-errors",
        dest="curate_annotation_errors",
        help=(
            "Apply activity curation based on pChEMBL values diverging in exactly 3.0, "
            "indicating possible annotation errors (e.g.: uM to nM and vice versa). Defaults to True."
        ),
        action="store_true",
    )
    parser.add_argument(
        "-rel",
        "--standard-relation",
        nargs="*",
        default=DEFAULTS["standard_relation"],
        help=("Filter only bioactivities with the specified relation. Not implemented yet."),
        type=str,
    )
    parser.add_argument(
        "-at",
        "--assay-types",
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
        choices=["trace", "debug", "info", "warning", "error", "critical"],
        help="Set the logging level. Defaults to info",
    )
    parser.add_argument(
        "-cr",
        "--chembl-release",
        type=int,
        help="chembl_release: specify latest ChEMBL release to extract data from (e.g., 28). Defaults to None.",
        default=DEFAULTS["chembl_release"],
    )
    parser.add_argument(
        "-mcols",
        "--metadata-columns",
        nargs="*",
        default=DEFAULTS["metadata_columns"],
        help="Extra metadata columns to keep in the final dataframe, aggregated by ';'. Defaults to [].",
        type=str,
    )
    parser.add_argument(
        "-idcols",
        "--id-columns",
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
        "-skip_agg",
        "--skip-not-aggregated",
        action="store_true",
        help="Skips saving the data before aggregation of same-molecule datapoint takes place.",
    )
    parser.add_argument(
        "-mutagg",
        "--aggregate-mutants",
        action="store_true",
        help=(
            "Aggregate data on targets regardless of the mutation flag in ChEMBL, treating "
            "mutants as the same target. Regardless of the configuration, mutation data will be "
            "stored under `mutation`. Default is False."
        ),
    )
    parser.add_argument(
        "-rec",
        "--skip-recipe",
        dest="skip_recipe",
        help=(
            "Skip saving the json file with the parameters used to fetch the data. "
            "The file will be saved with the asme output path, but with the `_recipe.json` suffix. "
            "Defaults to True."
        ),
        action="store_true",
    )
    parser.add_argument(
        "-back",
        "--chembl-backend",
        dest="chembl_backend",
        help=(
            "Choose the backend to use for interacting with the ChEMBL database. Options are 'webresource' "
            "and 'downloader'. If 'downloader' is used, the sql data will be downloaded using chembl_downloader. "
            "Using webresource will fetch the data using the ChEMBL webresource client instead, not requiring "
            "the download of the SQL database, but will be slower due to API requests. Defaults to 'downloader'."
        ),
        choices=["webresource", "downloader"],
        default=DEFAULTS["chembl_backend"],
        type=str,
    )
    parser.add_argument(
        "-v",
        "--chembl-version",
        dest="chembl_version",
        help=(
            "Not to be confused for the --chembl-release argument. Version of the ChEMBL database to be "
            "used by chembl_downloader. Only used when `backend=='downloader'`. If left as default, the "
            "latest version available will be used. Defaults to None."
        ),
        default=DEFAULTS["chembl_version"],
        type=str,
    )
    parser.add_argument(
        "-sd",
        "--save-dropped",
        dest="save_dropped",
        action="store_true",
        help=(
            "Save a separate dataframe containing rows that were flagged for dropping, "
            "along with the reasons for flagging. Default is False."
        ),
    )
    parser.add_argument(
        "-reqdoc",
        "--require-doc-date",
        dest="require_doc_date",
        action="store_true",
        help=("Filter out bioactivities that do not have a document date. Default is False."),
    )
    parser.add_argument(
        "-maxas",
        "--max-assay-size",
        dest="max_assay_size",
        type=int,
        default=DEFAULTS["max_assay_size"],
        help=(
            "Maximum number of compounds in an assay. Assays exceeding this size will have their "
            "activities flagged for removal. If left as default (None), this filter won't be applied."
        ),
    )
    parser.add_argument(
        "-minas",
        "--min-assay-size",
        dest="min_assay_size",
        type=int,
        default=DEFAULTS["min_assay_size"],
        help=(
            "Minimum number of compounds in an assay. Assays exceeding this size will have their "
            "activities flagged for removal. If left as default (None), this filter won't be applied."
        ),
    )
    parser.add_argument(
        "-maxm",
        "--max-assay-match",
        dest="max_assay_match",
        action="store_true",
        help=(
            "Perform assay metadata matching. If True, assay metadata columns will be added to "
            "`id_columns`, preventing compounds not matching the assay metadata from being "
            "aggregated. Assay metadata columns that define this matching are: "
            f"{', '.join(DEFAULT_ASSAY_MATCH_FIELDS)}. Default is False"
        ),
    )
    parser.add_argument(
        "-maso",
        "--min-assay-overlap",
        dest="min_assay_overlap",
        type=int,
        default=DEFAULTS["min_assay_overlap"],
        help=(
            "Minimum number of overlapping compounds between two assays for the same target "
            "for their activities to be considered. Activities from assay pairs not meeting "
            "this overlap will be flagged for removal. If left as default (None), "
            "this filter won't be applied."
        ),
    )
    parser.add_argument(
        "-smr",
        "--strict-mutant-removal",
        dest="strict_mutant_removal",
        action="store_true",
        help=(
            "If True, assays with 'mutant', 'mutation', or 'variant' in their "
            "description will be flagged for removal. Default is False."
        ),
    )
    parser.add_argument(
        "-kfd",
        "--keep-flagged-data",
        dest="keep_flagged_data",
        action="store_true",
        help=(
            "If True, data points flagged for dropping will be retained in the main dataset. "
            "The `data_dropping_comment` column will still be populated. "
            "A warning will be logged. Default is False."
        ),
    )
    parser.add_argument(
        "-cpd-eq",
        "--compound-equality",
        dest="compound_equality",
        choices=["mixed_fp", "connectivity"],
        default=DEFAULTS["compound_equality"],
        type=str,
        help=(
            "Method for compound equality determination. In case of 'mixed_fp', a mixed fingerprint "
            "composed of ECFP4 and RDKit fingerprints will be used for determining compound equality. "
            "Defaults to 'connectivity'."
        ),
    )

    return parser.parse_args()


def main(args: argparse.Namespace) -> None:

    setup_logger(level=args.log_level.upper())

    if args.command:
        if args.command == "download":
            check_and_download_chembl_db(prefix=args.prefix, version=args.version)
        elif args.command == "explore":
            explorer_main(
                version=args.version,
                list_tables=args.list_tables,
                table=args.table,
                search_column=args.search_column,
                query=args.query,
            )
        sys.exit(0)

    if args.standard_relation != ["="]:
        raise NotImplementedError("Fetching data using different relation types isn't implemented yet.")

    output_path = Path(args.output_path)
    if not output_path.parent.exists():
        output_path.mkdir()

    if args.chirality and not args.drop_unassigned_chiral:
        logger.warning(
            "Consider passing the `--drop_unassigned_chiral` flag when using the `--chirality` flag. "
            "For more information on why this could be problematic, see the link:\n"
            "https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00934-w#:~:text="
            "Some%20duplicates%20were,for%20kinetic%20solubility."
        )

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
        chembl_release=args.chembl_release,
        save_not_aggregated=(not args.skip_not_aggregated),
        save_dropped=args.save_dropped,
        drop_unassigned_chiral=args.drop_unassigned_chiral,
        curate_annotation_errors=args.curate_annotation_errors,
        version=args.chembl_version,
        backend=args.chembl_backend,
        require_doc_date=args.require_doc_date,
        min_assay_size=args.min_assay_size,
        max_assay_size=args.max_assay_size,
        min_assay_overlap=args.min_assay_overlap,
        strict_mutant_removal=args.strict_mutant_removal,
        keep_flagged_data=args.keep_flagged_data,
    )

    df = aggregate_data(
        df=df,
        chirality=args.chirality,
        extra_multival_cols=args.metadata_columns,
        extra_id_cols=args.id_columns,
        aggregate_mutants=args.aggregate_mutants,
        max_assay_match=args.max_assay_match,
        output_path=args.output_path,
        compound_equality=args.compound_equality,
    )

    # Save the recipe
    if not args.skip_recipe:
        output_name = output_path.stem
        recipe_path = output_path.parent / f"{output_name}_recipe.json"

        configs = vars(args)
        command_vals = []
        for k, v in configs.items():

            save_k = k.replace("_", "-")  # same format as the command line

            if k in ["chembl_release", "chembl_version"]:
                if v is not None:
                    command_vals.append((f"--{save_k} {v}"))
                else:  # safe to assume latest version; if None, chembl_downloader gets latest
                    command_vals.append(f"--{save_k} {latest()}")

            elif isinstance(v, bool):
                if DEFAULTS[k] is not v:
                    command_vals.append((f"--{save_k}" if k in STORE_TRUE_ARGS else f"--{k} {v}"))
            elif isinstance(v, list):
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{save_k} {' '.join([str(i) for i in v])}")
            elif isinstance(v, str):
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{save_k} {v}")
            elif isinstance(v, int):  # Added for max_assay_size and min_assay_overlap
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{save_k} {v}")

        command = "getchembl " + " ".join(command_vals)
        configs.update({"CompoundMapper version": __version__})
        configs.update({"command": command})

        with open(recipe_path, "w") as f:
            json.dump(configs, f, indent=2)

    return df


def main_exe() -> None:
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
