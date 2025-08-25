<div align="center">
  <img src="logo.svg" alt="" width=240>
  <p><strong>The ChEMBL data curator that flags issues instead of silently dropping them.</strong></p>

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Code style: black](https://img.shields.io/badge/code%20style-black-black?style=flat-square)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat-square&labelColor=ef8336)](https://pycqa.github.io/isort/)
[![License: MIT](https://img.shields.io/badge/License-MIT-purple?style=flat-square)](https://opensource.org/licenses/MIT)

</div>

> Inspired in the Portuguese word "*capricho*" [🔊](https://ipa-reader.com/?text=ka%CB%88p%C9%BEi.%CA%83u&voice=Ricardo). Doing someting *with capricho* means doing it *meticulously*, *with care* and *attention to detail*.

CAPRICHO (**C**hEMBL **A**ggregation **P**ackage with **R**obust **I**nspection and **C**uration **H**andling **O**ptions) is a Python package that streamlines fetching, curating, and aggregating ChEMBL data into a machine learning-ready format for drug discovery in a flexible and reproducible manner. Instead of making opiniated decisions on the source data, CAPRICHO curates it based on several quality control filters that can be chosen by the user. Its guiding principle is to never silently drop data. Entries that don't meet the criteria are marked, allowing the user to analyze how each curation step affects the comparability of assay readouts for the same compound.

## 🎯 Goals

The development of CAPRICHO is guided by two core principles:
- **Transparency Above All**: Data curation should never be a black box. Removed data points should be saved to be scrutinized by the user and the original data should be always preserved to ensure data integrity.
- **Flexibility by Design**: Every modeling project is unique. The tool must support flexible data collection and aggregation, allowing the incorporation of any ChEMBL metadata column to be incorporated into same-compound bioactivity values.

## ✨ Features:

- Data retrieval by any ChEMBL identifier (molecule IDs, target IDs, assay IDs, or document IDs)
- Automated pChEMBL (pXC50) value calculation for bioactivities if not provided through ChEMBL
- Customizable filtering options
- Configurable data aggregation options
- Save a fetching and processing recipe for reproducibility
- Command-line interface for easy use

## ⚙️ Installation

The most recent release can be installed from PyPI with uv:
```shell
uv pip install chembl_downloader
```

or with pip:
```shell
python -m pip install chembl_downloader
```

Alternatively, install directly from the GitHub repository with uv using the command:
```shell
uv pip install git+https://github.com/David-Araripe/Capricho.git
```
or with pip
```shell
python -m pip install git+https://github.com/David-Araripe/Capricho.git
```

## 🚀  Getting started

CAPRICHO provides a command-line interface with three main commands:
- [download](#download)
- [explore](#explore)
- [get](#get)

### Download

This command downloads the ChEMBL SQL database using `chembl_downloader`.

```bash
capricho download [OPTIONS]
```

**Options:**

| Option | Description | Default |
|---|---|---|
| `--version`, `-v` | ChEMBL version to download. | latest |
| `--prefix`, `-p` | Custom pystow storage path. | `~/.data/chembl/` |

### explore

Explore the downloaded ChEMBL SQL database.

```bash
capricho explore [OPTIONS]
```

**Options:**

| Option | Description |
|---|---|
| `--version`, `-v` | ChEMBL version to use. Defaults to the latest. |
| `--list-tables`, `-list` | List all tables within the SQL database and exit. |
| `--table`, `-t` | Explore a specific table. |
| `--search-column`, `-search` | Search for tables containing a column name pattern. |
| `--query`, `-q` | Run a custom SQL query. |


## capricho get

Filter, download, and process bioactivity data from ChEMBL. This is the main command of the package.

```bash
capricho get [OPTIONS]
```

A comprehensive list of options is provided below.

## Quickstart Examples

**1. Basic Target Query**

Get all bioactivity data for a specific target (e.g., EGFR) and save it to `egfr_data.csv`.

```bash
capricho get --target-ids CHEMBL203 --output-path egfr_data.csv
```

**2. Advanced Filtering**

Get data for a list of targets, but only with high confidence scores, and aggregate mutants.

```bash
capricho get --target-ids CHEMBL203,CHEMBL204 --confidence-scores 8,9 --aggregate-mutants --output-path advanced_query.csv
```

**3. Download a specific ChEMBL version**

Download ChEMBL version 33.

```bash
capricho download --version 33
```

## ⚙️ `get` Command Configuration Options

### Input IDs

| Option | Description | Default |
|---|---|---|
| `-mids`, `--molecule-ids` | ChEMBL molecule IDs, comma-separated. | `[]` |
| `-tids`, `--target-ids` | ChEMBL target IDs, comma-separated. | `[]` |
| `-asids`, `--assay-ids` | ChEMBL assay IDs, comma-separated. | `[]` |
| `-dids`, `--document-ids` | ChEMBL document IDs, comma-separated. | `[]` |

### Filtering Options

| Option | Description | Default |
|---|---|---|
| `-c`, `--confidence-scores` | Confidence scores to filter, comma-separated. | `[7, 8, 9]` |
| `-biotype`, `--bioactivity-type` | Bioactivity types to filter, comma-separated. | `['Potency', 'Kd', 'Ki', 'IC50', 'AC50', 'EC50']` |
| `-rel`, `--standard-relation` | Filter by standard relation, comma-separated. | `['=']` |
| `-at`, `--assay-types` | Assay types (B, F, A, T, P), comma-separated. | `['B', 'F']` |
| `-cr`, `--chembl-release` | Only fetch data reported **up to** a certain ChEMBL release. | `None` |
| `-reqdoc`, `--require-doc-date` | Filter out bioactivities without a document date. | `False` |
| `-maxas`, `--max-assay-size` | Max number of compounds in an assay. | `None` |
| `-minas`, `--min-assay-size` | Min number of compounds in an assay. | `None` |
| `-maso`, `--min-assay-overlap` | Min overlapping compounds between assays. | `0` |

### Processing & Aggregation Options

| Option | Description | Default |
|---|---|---|
| `-calc`, `--calculate-pchembl` | Calculate pChEMBL values if not reported. | `False` |
| `-chiral`, `--chirality` | Consider chirality during fingerprint calculation. | `False` |
| `-duchi`, `--drop-unassigned-chiral` | Drop entries with unassigned chiral centers. | `False` |
| `-cure`, `--curate-annotation-errors` | Apply curation for pChEMBL annotation errors. | `False` |
| `-mutagg`, `--aggregate-mutants` | Aggregate data on targets regardless of mutation. | `False` |
| `-maxm`, `--max-assay-match` | Perform strict assay metadata matching. | `False` |
| `-smr`, `--strict-mutant-removal` | Flag assays with mutant-related keywords for removal. | `False` |
| `-cpd-eq`, `--compound-equality` | Method for compound equality determination. | `connectivity` |
| `-mcols`, `--metadata-columns` | Extra metadata columns to keep, comma-separated. | `[]` |
| `-idcols`, `--id-columns` | Extra ID columns for aggregation, comma-separated. | `[]` |

### Output & Backend Options

| Option | Description | Default |
|---|---|---|
| `-o`, `--output-path` | Path to save the output files. | `chembl_data.csv` |
| `-skip-agg`, `--skip-not-aggregated` | Skip saving pre-aggregation data. | `False` |
| `-rec`, `--skip-recipe` | Skip saving the JSON recipe file. | `False` |
| `-back`, `--chembl-backend` | Backend to use for ChEMBL interaction. | `downloader` |
| `-v`, `--chembl-version` | ChEMBL version used by `chembl_downloader`. | `None` |

## Key Concepts

### Compound Equality

The `--compound-equality` option determines how CAPRICHO decides if two compound entries are the same.
- `connectivity`: (Default) Compounds are considered the same if they have the same chemical connectivity (i.e., ignoring stereochemistry).
- `mixed_fp`: A more stringent method using a combination of ECFP4 and RDKit fingerprints to determine similarity.

### ChEMBL Backends

The `--chembl-backend` option lets you choose how to fetch data:
- `downloader`: (Default) Uses a local SQL database downloaded via `chembl_downloader`. This is faster for large or repeated queries.
- `webresource`: Queries the live ChEMBL web API. This is useful for smaller, one-off queries without needing to download the entire database.

### Reproducibility with `recipe.json`

By default, CAPRICHO saves a `_recipe.json` file alongside your output data. This file contains the exact command and all the parameters you used to generate the data, ensuring your workflow is fully reproducible.

## Output Format

The fetcher returns a pandas DataFrame with the following key columns:
- `molecule_chembl_id`: ChEMBL ID for the compound
- `target_chembl_id`: ChEMBL ID for the target
- `standard_value`: Bioactivity measurement
- `standard_units`: Units of measurement
- `pchembl_value`: Calculated or reported pChEMBL value
- Additional metadata columns as specified

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.