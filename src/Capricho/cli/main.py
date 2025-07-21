"""Module containing the command line interface to get data from ChEMBL"""

import json
from enum import Enum
from pathlib import Path
from typing import List, Optional

import numpy as np
import typer
from typing_extensions import Annotated

from chembl_downloader import latest

from .. import __version__
from ..chembl.api.downloader import check_and_download_chembl_db
from ..chembl.api.sql_explorer import explorer_main
from ..core.default_fields import DEFAULT_ASSAY_MATCH_FIELDS
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
    "drop_unassigned_chiral": False,
    "curate_annotation_errors": False,
    "standard_relation": ["="],
    "assay_types": ["B", "F"],
    "log_level": "info",
    "chembl_release": None,
    "metadata_columns": [],
    "id_columns": [],
    "skip_not_aggregated": False,
    "aggregate_mutants": False,
    "skip_recipe": False,
    "chembl_backend": "downloader",
    "chembl_version": None,
    "require_doc_date": False,
    "max_assay_size": None,
    "min_assay_size": None,
    "max_assay_match": False,
    "min_assay_overlap": 0,
    "strict_mutant_removal": False,
    "compound_equality": "connectivity",
}

STORE_TRUE_ARGS = [
    "calculate_pchembl",
    "chirality",
    "drop_unassigned_chiral",
    "curate_annotation_errors",
    "skip_not_aggregated",
    "aggregate_mutants",
    "skip_recipe",
    "require_doc_date",
    "max_assay_match",
    "strict_mutant_removal",
]


app = typer.Typer(
    name="CAPRICHO",
    help="A ChEMBL data curator that flags questionable entries instead of silently dropping them.",
    no_args_is_help=True,
    rich_markup_mode="markdown",
    context_settings={"help_option_names": ["--help", "-h"], "max_content_width": 88},
)


def csv_string(value: Optional[str]) -> List[str]:
    """Parses a comma-separated string into a list of strings."""
    logger.debug(f"Parsing CSV string type= {type(value)}, value= {value}")
    if isinstance(value, list):
        return value
    elif value is None:
        raise ValueError("Value cannot be None")
    return [item.strip() for item in value.split(",")]


def csv_intergers(value: Optional[str]) -> List[int]:
    """Parses a comma-separated string into a list of integers."""
    logger.debug(f"Parsing CSV int type={type(value)}, value={value}")
    if isinstance(value, list):
        return [int(item) for item in value if isinstance(item, int)]
    elif isinstance(value, int):
        return [value]
    elif value is None:
        raise ValueError("Value cannot be None")
    return [int(item.strip()) for item in value.split(",") if item.strip().isdigit()]


class LogLevel(str, Enum):
    trace = "trace"
    debug = "debug"
    info = "info"
    warning = "warning"
    error = "error"
    critical = "critical"


class ChemblBackend(str, Enum):
    downloader = "downloader"
    webresource = "webresource"


class CompoundEquality(str, Enum):
    mixed_fp = "mixed_fp"
    connectivity = "connectivity"


@app.callback()
def main(
    ctx: typer.Context,  # noqa: F821
    log_level: Annotated[
        LogLevel,
        typer.Option(
            "-log",
            "--log-level",
            help="Set the logging level.",
            case_sensitive=False,
        ),
    ] = DEFAULTS.get("log_level", "info"),
):
    """
    Manage the ChEMBL data-wrangling workflow.
    """
    ctx.obj = {}
    logger.info(f"Log level set to {log_level.value.upper()}")
    setup_logger(level=log_level.upper())


@app.command()
def download(
    ctx: typer.Context,  # noqa: F821
    version: Annotated[
        Optional[int],
        typer.Option("--version", "-v", help="ChEMBL version to download. Defaults to the latest."),
    ] = None,
    prefix: Annotated[
        Optional[Path],
        typer.Option("--prefix", "-p", help="Custom pystow storage path. Defaults to ~/.data/chembl/."),
    ] = None,
):
    """Download ChEMBL SQL database using chembl_downloader."""
    logger.info(f"Starting ChEMBL download command for version: {version or 'latest'}")
    check_and_download_chembl_db(prefix=str(prefix) if prefix else None, version=version)
    raise typer.Exit()


@app.command()
def explore(
    ctx: typer.Context,  # noqa: F821
    version: Annotated[
        Optional[int], typer.Option("--version", "-v", help="ChEMBL version to use. Defaults to the latest.")
    ] = None,
    list_tables: Annotated[
        bool, typer.Option("--list-tables", "-list", help="List all tables within the SQL database and exit.")
    ] = False,
    table: Annotated[Optional[str], typer.Option("--table", "-t", help="Explore a specific table.")] = None,
    search_column: Annotated[
        Optional[str],
        typer.Option(
            "--search-column", "-search", help="Search for tables containing a column name pattern."
        ),
    ] = None,
    query: Annotated[Optional[str], typer.Option("--query", "-q", help="Run a custom SQL query.")] = None,
):
    """
    Explore the downloaded ChEMBL SQL database.

    For a visual inspection of the latest ChEMBL schema, see: https://www.ebi.ac.uk/chembl/db_schema
    """
    logger.info("Starting ChEMBL explore command.")
    explorer_main(
        version=version,
        list_tables=list_tables,
        table=table,
        search_column=search_column,
        query=query,
    )
    raise typer.Exit()


@app.command(name="get", no_args_is_help=True)
def get_data(
    ctx: typer.Context,
    # --- Input ID Arguments ---
    molecule_ids: Annotated[
        Optional[str],
        typer.Option(
            "-mids",
            "--molecule-ids",
            parser=csv_string,
            help="ChEMBL molecule IDs, comma-separated.",
            show_default=True,
            metavar="ID,ID,...",
        ),
    ] = DEFAULTS["molecule_ids"],
    target_ids: Annotated[
        Optional[str],
        typer.Option(
            "-tids",
            "--target-ids",
            parser=csv_string,
            help="ChEMBL target IDs, comma-separated.",
            show_default=True,
            metavar="ID,ID,...",
        ),
    ] = DEFAULTS["target_ids"],
    assay_ids: Annotated[
        Optional[str],
        typer.Option(
            "-asids",
            "--assay-ids",
            parser=csv_string,
            help="ChEMBL assay IDs, comma-separated.",
            show_default=True,
            metavar="ID,ID,...",
        ),
    ] = DEFAULTS["assay_ids"],
    document_ids: Annotated[
        Optional[str],
        typer.Option(
            "-dids",
            "--document-ids",
            parser=csv_string,
            help="ChEMBL document IDs, comma-separated.",
            show_default=True,
            metavar="ID,ID,...",
        ),
    ] = DEFAULTS["document_ids"],
    # --- Configuration Arguments ---
    output_path: Annotated[
        Path,
        typer.Option(
            "-o",
            "--output-path",
            help="Path to save the output files.",
        ),
    ] = DEFAULTS["output_path"],
    confidence_scores: Annotated[
        str,
        typer.Option(
            "-c",
            "--confidence-scores",
            help="Confidence scores to filter, comma-separated.",
            parser=csv_intergers,
            metavar="SCORES[1-9]",
            show_default=True,
        ),
    ] = DEFAULTS["confidence_scores"],
    bioactivity_type: Annotated[
        str,
        typer.Option(
            "-biotype",
            "--bioactivity-type",
            parser=csv_string,
            help="Bioactivity types to filter, comma-separated.",
            metavar="Ki,Kd,...",
        ),
    ] = DEFAULTS["bioactivity_type"],
    standard_relation: Annotated[
        str,
        typer.Option(
            "-rel",
            "--standard-relation",
            parser=csv_string,
            help="Filter by standard relation, comma-separated.",
            metavar="=,>,<",
        ),
    ] = DEFAULTS["standard_relation"],
    assay_types: Annotated[
        str,
        typer.Option(
            "-at",
            "--assay-types",
            parser=csv_string,
            help="Assay types (B, F, A, T, P), comma-separated.",
            metavar="B,F",
        ),
    ] = DEFAULTS["assay_types"],
    chembl_release: Annotated[
        Optional[int],
        typer.Option(
            "-cr",
            "--chembl-release",
            help="Only fetch data reported **up to** a certain ChEMBL release (e.g., 34).",
            metavar="int",
        ),
    ] = DEFAULTS["chembl_release"],
    chembl_version: Annotated[
        Optional[str],
        typer.Option(
            "-v",
            "--chembl-version",
            help="ChEMBL version used by _chembl_downloader_.",
            metavar="str",
        ),
    ] = DEFAULTS["chembl_version"],
    chembl_backend: Annotated[
        ChemblBackend,
        typer.Option(
            "-back",
            "--chembl-backend",
            help="Backend to use for ChEMBL interaction.",
        ),
    ] = DEFAULTS["chembl_backend"],
    compound_equality: Annotated[
        CompoundEquality,
        typer.Option(
            "-cpd-eq",
            "--compound-equality",
            help="Method for compound equality determination. mixed_fp uses combined ECFP4 and RDKit fingerprints.",
        ),
    ] = DEFAULTS["compound_equality"],
    # --- Metadata & Aggregation ---
    metadata_columns: Annotated[
        List[str],
        typer.Option(
            "-mcols",
            "--metadata-columns",
            parser=csv_string,
            help="Extra metadata columns to keep, comma-separated.",
            show_default=False,
            metavar="col1,col2,...",
        ),
    ] = DEFAULTS["metadata_columns"],
    id_columns: Annotated[
        List[str],
        typer.Option(
            "-idcols",
            "--id-columns",
            parser=csv_string,
            help="Extra ID columns for aggregation, comma-separated. E.g.: 'assay_chembl_id'",
            show_default=False,
            metavar="col1,col2,...",
        ),
    ] = DEFAULTS["id_columns"],
    # --- Numeric Filter Arguments ---
    max_assay_size: Annotated[
        Optional[int],
        typer.Option(
            "-maxas",
            "--max-assay-size",
            help="Max number of compounds in an assay.",
            metavar="int",
        ),
    ] = DEFAULTS["max_assay_size"],
    min_assay_size: Annotated[
        Optional[int],
        typer.Option(
            "-minas",
            "--min-assay-size",
            help="Min number of compounds in an assay.",
            metavar="int",
        ),
    ] = DEFAULTS["min_assay_size"],
    min_assay_overlap: Annotated[
        int,
        typer.Option(
            "-maso",
            "--min-assay-overlap",
            help="Min overlapping compounds between assays.",
            metavar="int",
        ),
    ] = DEFAULTS["min_assay_overlap"],
    # --- Boolean Flags ---
    calculate_pchembl: Annotated[
        bool,
        typer.Option(
            "-calc/-no-calc",
            "--calculate-pchembl/--no-calculate-pchembl",
            help="Calculate pChEMBL values if not reported.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["calculate_pchembl"],
    chirality: Annotated[
        bool,
        typer.Option(
            "-chiral/-no-chiral",
            "--chirality/--no-chirality",
            help="Consider chirality during fingerprint calculation.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["chirality"],
    drop_unassigned_chiral: Annotated[
        bool,
        typer.Option(
            "-duchi/-dont-duchi",
            "--drop-unassigned-chiral/--dont-drop-unassigned-chiral",
            help="Drop entries with unassigned chiral centers.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["drop_unassigned_chiral"],
    curate_annotation_errors: Annotated[
        bool,
        typer.Option(
            "-cure/-dont-cure",
            "--curate-annotation-errors/--dont-curate-annotation-errors",
            help="Apply curation for pChEMBL annotation errors.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["curate_annotation_errors"],
    skip_not_aggregated: Annotated[
        bool,
        typer.Option(
            "-skip-agg/-dont-skip-agg",
            "--skip-not-aggregated/--dont-skip-not-aggregated",
            help="Skip saving pre-aggregation data.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["skip_not_aggregated"],
    aggregate_mutants: Annotated[
        bool,
        typer.Option(
            "-mutagg/-dont-mutagg",
            "--aggregate-mutants/--dont-aggregate-mutants",
            help="Aggregate data on targets regardless of mutation.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["aggregate_mutants"],
    skip_recipe: Annotated[
        bool,
        typer.Option(
            "-rec/-dont-rec",
            "--skip-recipe/--dont-skip-recipe",
            help="Skip saving the JSON recipe file.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["skip_recipe"],
    require_doc_date: Annotated[
        bool,
        typer.Option(
            "-reqdoc/-dont-reqdoc",
            "--require-doc-date/--dont-require-doc-date",
            help="Filter out bioactivities without a document date.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["require_doc_date"],
    max_assay_match: Annotated[
        bool,
        typer.Option(
            "-maxm/-no-maxm",
            "--max-assay-match/--no-max-assay-match",
            help=f"Perform strict assay metadata matching based on: {', '.join(DEFAULT_ASSAY_MATCH_FIELDS)}.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["max_assay_match"],
    strict_mutant_removal: Annotated[
        bool,
        typer.Option(
            "-smr/-no-smr",
            "--strict-mutant-removal/--dont-strict-mutant-removal",
            help="Flag assays with mutant-related keywords in assay_description for removal.",
            is_flag=True,
            metavar="bool",
        ),
    ] = DEFAULTS["strict_mutant_removal"],
):
    """
    Filter, download, and process bioactivity data from ChEMBL.
    """
    if standard_relation != ["="]:
        logger.error("Fetching data using different relation types isn't implemented yet.")
        raise typer.Exit(code=1)

    if not output_path.parent.exists():
        output_path.mkdir()
    if output_path.suffix == "":
        output_path = output_path.with_suffix(".csv")

    valid_suffixes = [".csv", ".tsv", ".parquet"]
    if not np.intersect1d(valid_suffixes, output_path.suffixes).shape[0] > 0:
        logger.error(
            "Output file must have a valid suffix: .csv, .tsv, or .parquet. "
            f"Provided suffix: {output_path.suffixes}"
        )
        raise typer.Exit(code=1)

    if chirality and not drop_unassigned_chiral:
        logger.warning(
            "Consider passing the `--drop-unassigned-chiral` flag when using the `--chirality` flag. "
            "For more information on why this could be problematic, see the link:\n"
            "https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00934-w#:~:text="
            "Some%20duplicates%20were,for%20kinetic%20solubility."
        )

    df = get_standardize_and_clean_workflow(
        molecule_ids=molecule_ids or [],
        target_ids=target_ids or [],
        assay_ids=assay_ids or [],
        document_ids=document_ids or [],
        chirality=chirality,
        calculate_pchembl=calculate_pchembl,
        output_path=output_path,
        confidence_scores=confidence_scores,
        bioactivity_type=bioactivity_type,
        standard_relation=standard_relation,
        assay_types=assay_types,
        chembl_release=chembl_release,
        save_not_aggregated=(not skip_not_aggregated),
        drop_unassigned_chiral=drop_unassigned_chiral,
        curate_annotation_errors=curate_annotation_errors,
        version=chembl_version,
        backend=chembl_backend.value,
        require_doc_date=require_doc_date,
        min_assay_size=min_assay_size,
        max_assay_size=max_assay_size,
        min_assay_overlap=min_assay_overlap,
        strict_mutant_removal=strict_mutant_removal,
    )

    df = aggregate_data(
        df=df,
        chirality=chirality,
        extra_multival_cols=metadata_columns,
        extra_id_cols=id_columns,
        aggregate_mutants=aggregate_mutants,
        max_assay_match=max_assay_match,
        output_path=output_path,
        compound_equality=compound_equality.value,
    )

    if not skip_recipe:
        output_name = output_path.stem.split(".")[0]
        recipe_path = output_path.parent / f"{output_name}_recipe.json"

        command_vals = []
        configs = {k: v for k, v in ctx.params.items() if v is not None}
        for k, v in configs.items():
            save_k = k.replace("_", "-")  # same format as the command line
            if k in ["chembl_release", "chembl_version"]:
                if v is not None:
                    command_vals.append((f"--{save_k} {v}"))
                else:  # safe to assume latest version; if None, chembl_downloader gets latest
                    command_vals.append(f"--{save_k} {latest()}")
            elif isinstance(v, Path) or isinstance(v, Enum):
                configs[k] = str(v)  # Convert Path and Enum to string for JSON serialization
                if DEFAULTS[k] is not v:
                    command_vals.append(f"--{save_k} {str(v)}")
            elif isinstance(v, bool):
                if DEFAULTS[k] is not v:
                    command_vals.append((f"--{save_k}" if k in STORE_TRUE_ARGS else f"--{k} {v}"))
            elif isinstance(v, list):
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{save_k} {','.join([str(i) for i in v])}")
            elif isinstance(v, str):
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{save_k} {v}")
            elif isinstance(v, int):
                if DEFAULTS[k] != v:
                    command_vals.append(f"--{save_k} {v}")
        command = "capricho get " + " ".join(command_vals)
        configs = {k: configs[k] for k in DEFAULTS if k in configs}
        configs = {"command": command, "capricho version": __version__, **configs}

        with open(recipe_path, "w") as f:
            json.dump(configs, f, indent=2)

        logger.info(f"Recipe saved to {recipe_path}")

    logger.info(f"Successfully processed and saved data to {output_path}")
    return df


if __name__ == "__main__":
    app()
