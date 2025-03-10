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
| `chembl_version` | Specific ChEMBL version to use | `None` |
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
| `aggregate_mutants` | Aggregate data on targets regardless of their variant sequence, treating mutants as the same target. Mutation data is still stored under `variant_sequence` in ChEMBL | `False` |
| `skip_not_aggregated` | Skip saving the raw data before any aggregation of repeated molecules is applied. Not-aggregated data is saved with the `_not_aggregated.csv` suffix | `False` |
| `calculate_pchembl` | Calculate pChEMBL (pXC50) values for bioactivities reported in nM, µM or uM when not available | `False` |
| `no_document_info` | Skip retrieving ChEMBL document info to reduce API calls. Passing this has the drawback that information such as `year` and `chembl_release` will be missing  | `False` |
| `drop_unassigned_chiral` | Drop ChEMBL compounds with one or more undefined stereocenters. Advised when passing the `-chiral` flag | `False` |

### Quality filters

Several quality filters have been previously mentioned in the literature and are supported in CompoundMapper, granting the user full control over the aggregated datsets' quality. We attribute the filters to the scientific papers that have used them:

- Quality filters from [Landrum & Riniker, 2024](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049).

| Filter Name | Explanation | Applying in CompoundMapper |
|-------------|-------------|---------|
| `activity_curation` | Remove pairs of measurements where pchembl values in two assays were either exactly the same or differed by 3.0 | Applied by default upon data curation |
| `duplicate_papers` | Remove measurements where both assays were published in the same document. Usually only occurs when there is a difference between the two assays | Passing parameters to the `--id_columns` parameter should prevent such undesired aggregation from taking place |
| `remove_mutants` | Remove any assays that have any of the keywords "mutant", "mutation", or "variant" keywords in the assay description. | In recent ChEMBL releases, such information is stored under the `variant_sequence` field, saved to the metadata. The conservative keyword search is also available using the `--conservative_remove_mutants` argument |
| `assay_metadata`| Remove pairs of assays with not-matching assay metadata such as *assay_type*, *assay_organism*, *assay_category*, *assay_tax_id*, *assay_strain*, *assay_tissue*, *assay_cell_type*, *assay_subcellular_fraction*, and *bao_format* | Aggregate only data points from matching assays by simply passing such values to the `--id_columns` argument |
| `sources_other_than_documents`| Remove data points that does not have an associated document date. E.g.: screening data sets or other contributed data sets. | Pass the argument `--assay_documents_only` |
| `assay_size`| remove assays that have `N > value` reported compounds. Used in the paper to focus on primary literature, avoiding review articles | Available through the argument `--assay_max_size` |
| `curation_confidence`| only query assays that have a certain confidence score (used only assays with confidence 9 in the paper) | available through the `--confidence_scores` argument |

<!-- Curation methods from https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049:

Activity curation: Pairs of measurements where the pchembl values in the two assays were either exactly the same or differed by 3.0 were removed. Given the very low probability of two separate experiments producing exactly the same results, the exact matches are most likely cases where values from a previous paper are copied into a new one; this was discussed in the earlier work by Kramer et al. (10) and spot-checked with a number of assay pairs here.
Duplicate papers: ✅
Remove mutants: ✅
Assay type: This curation step removes pairs of assays with different assay types. ✅
Assay metadata: ❌❌❌ (TODO): this curation step removes pairs of assays where any of the following assay metadata fields do not match: assay_type, assay_organism, assay_category, assay_tax_id, assay_strain, assay_tissue, assay_cell_type, assay_subcellular_fraction, and bao_format.
Sources other than documents: ❌❌❌ (TODO) this curation step removes any assay that is from a source that does not have an associated document date. 
Assay size ❌❌❌ (TODO): by default, any assays that include >100 compounds are removed.
Curation confidence: ✅
 -->

- Curation filters from [Tetko et al. 2024](https://doi.org/10.1186/s13321-024-00934-w).

| Filter Name | Explanation | Applying in CompoundMapper |
|-------------|-------------|---------|
| `remove_unassigned_chiral` | Authors found duplicated data points within a dataset due to compounds being present with defined and undefined stereoceters. | When modeling datasets with stereochemistry, avoid this problem by passing the `--drop_unassigned_chiral` argument |


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
