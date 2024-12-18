# CompoundMapper

A Python Package for efficiently curating bioactivity data from the ChEMBL database.

Simplify fetching, standardizing, and aggregating bioactivity data, outputting a machine learning-ready dataset for drug discovery that can be shared, reproduced, and updated at any ChEMBL release.

## Features:

- Flexible data retrieval by molecule IDs, target IDs, assay IDs, or document IDs
- Automated pChEMBL (pXC50) value calculation for bioactivities if not provided through ChEMBL
- Customizable filtering options (see below):
    - Confidence score filtering
    - Bioactivity type selection (Potency, Kd, Ki, IC50, AC50, EC50)
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
| `chembl_version` | Specific ChEMBL version to use | `None` |
| `standard_relation` | ChEMBL standard relations to use | [`=`] |

### 🔨 Processing Options

Once the data is retrieved, CompoundMapper executes the following processing steps:
- Merge all the retrieved data into a single DataFrame
- Standardize the SMILES strings using the [ChEMBL_Structure_Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline) package
- Drop entries with missing SMILES strings
- Remove salts and solvent molecules from the SMILES strings using a regex pattern defined in `CompoundMapper.core.smiles_utils.MIXTURE_REGEX`

The other operations are performed based on the following configuration options:

| Option | Description | Default |
|--------|-------------|---------|
| `id_columns` | Additional columns to use as identifiers during aggregation. For example, using `assay_chembl_id` will only aggregate data if both the compound and assay are identical | `[]` |
| `chirality` | Consider chirality when calculating molecular fingerprints. Important for differentiating between enantiomers during data aggregation | `False` |
| `aggregate_mutants` | Aggregate data on targets regardless of their variant sequence, treating mutants as the same target. Mutation data is still stored under `variant_sequence` in ChEMBL | `False` |
| `save_not_aggregated` | Save the raw data before performing any aggregation of repeated molecules | `False` |
| `calculate_pchembl` | Calculate pChEMBL (pXC50) values for bioactivities reported in nM, µM or uM when not available | `False` |
| `no_document_info` | Skip retrieving document information (like publication year) to reduce API calls. Passing this has the drawback that informations such as `year` and `chembl_release` will be missing  | `False` |


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
