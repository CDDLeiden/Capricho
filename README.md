# CompoundMapper

A Python Package for efficiently curating bioactivity data from the ChEMBL database.

Simplify fetching, standardizing, and aggregating bioactivity data, outputting a machine learning-ready dataset for drug discovery that can be shared, reproduced, and updated at any ChEMBL release.

## Features:

- Flexible data retrieval by molecule IDs, target IDs, assay IDs, or document IDs
- Automated pChEMBL (pXC50) value calculation for bioactivities if not provided through ChEMBL
- Customizable filtering options (see below):
    - Confidence score filtering
    - Bioactivity type selection (e.g.: Potency, Kd, Ki, IC50, AC50, EC50)
    - Assay type filtering (Functional, Binding, ADME, Toxicity, Physicochemical)
    - Standard relation filtering
- Configurable data aggregation options (see below)
- Support for multiple ChEMBL versions
- Recipe saving for reproducibility

Further, the package also methods for:
- Multithreaded compound structure curation through pubchempy
- Multithreaded compound ChEMBL similarity search

## ⚙️ Configuration Options

Retrieving data from ChEMBL is supported by any of the starting points:
- `molecule_ids`
- `target_ids`
- `assay_ids`
- `document_ids`

### 🔍 Filtering Options

Once you have the starting point, additional parameters can be passed to filter the data:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `confidence_scores` | Confidence scores for filtering bioactivities | [`7`, `8`, `9`] |
| `bioactivity_type` | Types of bioactivity to include | [`Potency`, `Kd`, `Ki`, `IC50`, `AC50`, `EC50`] |
| `assay_types` | Types of assays to include | [`B`, `F`] |
| `chembl_release` | Get only activities up to a certain release. `None` gets up to the latest reported activities. | `None` |
| `standard_relation` | ChEMBL standard relations to use | [`=`] |

### 🔨 Processing Options

Once the data is retrieved, CompoundMapper executes the following processing steps:
- Merge all the retrieved data into a single DataFrame
- Standardize the SMILES strings using the [ChEMBL_Structure_Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline) package
- Drop entries with missing SMILES strings
- Remove salts and solvent molecules from the SMILES strings using a regex pattern defined in `CompoundMapper.core.smiles_utils.MIXTURE_REGEX`
- Drop any remaining mixtures (e.g., multiple compounds in a single SMILES string separated by a dot)

The other operations performed to aggregate the fetched dataset and produce the ML-ready dataset are configurable based on the following options:

| Option | Description | Default |
|--------|-------------|---------|
| `id_columns` | Additional columns to use as identifiers during aggregation. For example, using `assay_chembl_id` will **only** aggregate readouts if the assay reporting them is the same | `[]` |
| `chirality` | Consider whether molecules have the same stereochemistry when aggregating repeated bioactivity measurements | `False` |
| `aggregate_mutants` | Aggregate data on targets regardless of their variant sequence, treating mutants as the same target. Mutation data is still stored under `mutation` in ChEMBL | `False` |
| `skip_not_aggregated` | Skip saving the raw data before any aggregation of repeated molecules is applied. Not-aggregated data is saved with the `_not_aggregated.csv` suffix | `False` |
| `calculate_pchembl` | Calculate pChEMBL (pXC50) values for bioactivities reported in nM, µM or uM when not available | `False` |
| `no_document_info` | Skip retrieving ChEMBL document info to reduce API calls. Passing this has the drawback that information such as `year` and `chembl_release` will be missing  | `False` |
| `drop_unassigned_chiral` | Drop ChEMBL compounds with one or more undefined stereocenters. Advised when passing the `-chiral` flag | `False` |

<!-- TODO: add the link to the `drop_unassigned_chiral` option: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00934-w#:~:text=Some%20duplicates%20were,for%20kinetic%20solubility. -->

### Output Options


## Output Format

The fetcher returns a pandas DataFrame with the following key columns:
- `molecule_chembl_id`: ChEMBL ID for the compound
- `target_chembl_id`: ChEMBL ID for the target
- `standard_value`: Bioactivity measurement
- `standard_units`: Units of measurement
- `pchembl_value`: Calculated or reported pChEMBL value
- Additional metadata columns as specified

## Installation:
Clone the repo and then you can install it in your envionment by:
```bash
python -m pip install -e .
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
