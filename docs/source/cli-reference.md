# CLI Reference

CAPRICHO provides four main commands: `download`, `explore`, `get`, and `binarize`. This section provides comprehensive documentation for all command-line options.

## capricho download

Downloads the ChEMBL SQL database using `chembl_downloader`.

```bash
capricho download [OPTIONS]
```

### Options

| Option | Description | Default |
|---|---|---|
| `--version`, `-v` | ChEMBL version to download | latest |
| `--prefix`, `-p` | Custom pystow storage path | `~/.data/chembl/` |

### Examples

```bash
# Download latest ChEMBL version
capricho download

# Download specific version
capricho download --version 33

# Use custom storage location (this will install version 25 it on ~/.data/old-chembl)
capricho download --version 25 --prefix data/old-chembl/
```

## capricho explore

Explore the downloaded ChEMBL SQL database.

```bash
capricho explore [OPTIONS]
```

### Options

| Option | Description |
|---|---|
| `--version`, `-v` | ChEMBL version to use (defaults to latest) |
| `--list-tables`, `-list` | List all tables within the SQL database and exit |
| `--table`, `-t` | Explore a specific table |
| `--search-column`, `-search` | Search for tables containing a column name pattern |
| `--query`, `-q` | Run a custom SQL query |

### Examples

```bash
# List all available tables on the relational database
capricho explore --list-tables

# Examine the activities table
capricho explore --table activities

# Find tables with 'pchembl' columns
capricho explore --search-column pchembl

# Check available standard_relation types in ChEMBL
capricho explore --query "SELECT standard_relation, COUNT(*) as count FROM activities WHERE standard_relation IS NOT NULL GROUP BY standard_relation ORDER BY count DESC"
```

### Useful Queries for Data Exploration

These queries can help you understand the data before running `capricho get`:

**Count activities by bioactivity type:**
```bash
capricho explore --query "SELECT standard_type, COUNT(*) as count FROM activities WHERE standard_type IS NOT NULL GROUP BY standard_type ORDER BY count DESC LIMIT 20"
```

**Count activities by assay type:**
```bash
capricho explore --query "SELECT assay_type, COUNT(*) as count FROM assays GROUP BY assay_type ORDER BY count DESC"
```

**Check confidence score distribution:**
```bash
capricho explore --query "SELECT confidence_score, COUNT(*) as count FROM assays GROUP BY confidence_score ORDER BY confidence_score DESC"
```

## capricho get

Filter, download, and process bioactivity data from ChEMBL. This is the main command of CAPRICHO.

```bash
capricho get [OPTIONS]
```

### Input ID Options

Specify which ChEMBL entities to retrieve data for:

| Option | Description | Default |
|---|---|---|
| `-mids`, `--molecule-ids` | ChEMBL molecule IDs, comma-separated | `[]` |
| `-tids`, `--target-ids` | ChEMBL target IDs, comma-separated | `[]` |
| `-asids`, `--assay-ids` | ChEMBL assay IDs, comma-separated | `[]` |
| `-dids`, `--document-ids` | ChEMBL document IDs, comma-separated | `[]` |

### Filtering Options

Control which bioactivity data to include:

| Option | Description | Default |
|---|---|---|
| `-c`, `--confidence-scores` | Confidence scores to filter, comma-separated | `[7, 8, 9]` |
| `-biotype`, `--bioactivity-type` | Bioactivity types to filter, comma-separated | `['Potency', 'Kd', 'Ki', 'IC50', 'AC50', 'EC50']` |
| `-rel`, `--standard-relation` | Filter by standard relation (`=`, `<`, `>`, `~`), comma-separated. **Note:** Including `<` or `>` requires `--calculate-pchembl`. See [Standard Relations](concepts.md). | `['=']` |
| `-at`, `--assay-types` | Assay types (B, F, A, T, P), comma-separated | `['B', 'F']` |
| `-cr`, `--chembl-release` | Only fetch data reported **up to** a certain ChEMBL release | `None` |
| `-reqdoc`, `--require-doc-date` | Filter out bioactivities without a document date | `False` |
| `-maxas`, `--max-assay-size` | Maximum number of compounds in an assay | `None` |
| `-minas`, `--min-assay-size` | Minimum number of compounds in an assay | `None` |
| `-maso`, `--min-assay-overlap` | Minimum overlapping compounds between assays | `0` |

#### Confidence Scores
ChEMBL assigns confidence scores from 0-9:
- **9**: Direct single protein target
- **8**: Direct protein complex/family target  
- **7**: Direct target with homologues
- **6**: Molecular target assigned
- **5**: Non-molecular target assigned
- **4**: Subcellular target assigned
- **3**: Cell-line target assigned
- **2**: Tissue target assigned
- **1**: Organism target assigned
- **0**: Unchecked data

#### Bioactivity Types
Common bioactivity measurements:
- **IC50**: Half maximal inhibitory concentration
- **EC50**: Half maximal effective concentration
- **Ki**: Inhibition constant
- **Kd**: Dissociation constant
- **AC50**: Half maximal activity concentration
- **Potency**: General potency measurement (see `assay_description` for more detail)

#### Assay Types
- **B**: Binding assay
- **F**: Functional assay
- **A**: ADMET assay
- **T**: Toxicity assay
- **P**: Physicochemical assay
- **U**: Unclassified

Further information on assay types and confidence scores can be found in the [ChEMBL documentation](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/chembl-data-questions).

### Processing & Aggregation Options

Control how data is processed and aggregated:

| Option | Description | Default |
|---|---|---|
| `-calc`, `--calculate-pchembl` | Calculate pChEMBL values if not reported. **Required when using censored data** (`--standard-relation` includes `<` or `>`). See [Standard Relations](concepts.md). | `False` |
| `-chiral`, `--chirality` | Consider chirality during fingerprint calculation | `False` |
| `-duchi`, `--drop-unassigned-chiral` | Drop entries with unassigned chiral centers | `False` |
| `-cure`, `--curate-annotation-errors` | Apply curation for pChEMBL annotation errors | `False` |
| `-mutagg`, `--aggregate-mutants` | Aggregate data on targets regardless of mutation | `False` |
| `-maxm`, `--max-assay-match` | Perform strict assay metadata matching | `False` |
| `-smr`, `--strict-mutant-removal` | Flag assays with mutant-related keywords for removal | `False` |
| `-cpd-eq`, `--compound-equality` | Method for compound equality determination | `connectivity` |
| `-mcols`, `--metadata-columns` | Extra metadata columns to keep, comma-separated | `[]` |
| `-idcols`, `--id-columns` | Extra ID columns for aggregation, comma-separated | `[]` |

#### Compound Equality Methods
- **connectivity**: (Default) Based on molecular connectivity, ignoring stereochemistry
- **mixed_fp**: Uses ECFP4 and RDKit fingerprints (each with 2048 bits) for similarity determination

#### Useful Metadata Columns
- `organism`: Source organism
- `tissue`: Tissue type
- `cell_type`: Cell line information
- `assay_description`: Detailed assay description
- `target_type`: Type of target (e.g., SINGLE PROTEIN, PROTEIN FAMILY)

### Output & Backend Options

Control output format and data source:

| Option | Description | Default |
|---|---|---|
| `-o`, `--output-path` | Path to save the output files | `chembl_data.csv` |
| `-skip-agg`, `--skip-not-aggregated` | Skip saving pre-aggregation data | `False` |
| `-rec`, `--skip-recipe` | Skip saving the JSON recipe file | `False` |
| `-back`, `--chembl-backend` | Backend to use for ChEMBL interaction | `downloader` |
| `-v`, `--chembl-version` | ChEMBL version used by `chembl_downloader` | `None` |

#### Backend Options
- **downloader**: (Default) Uses local SQL database, faster for large queries
- **webresource**: Uses ChEMBL web API, no local download required

## Complete Examples

### Basic Target Analysis
```bash
capricho get --target-ids CHEMBL203 --output-path egfr_analysis.csv
```

### High-Quality Multi-Target Dataset
```bash
capricho get \
  --target-ids CHEMBL203,CHEMBL204,CHEMBL279 \
  --confidence-scores 8,9 \
  --bioactivity-type IC50,Ki \
  --standard-relation "=" \
  --aggregate-mutants \
  --metadata-columns organism,tissue \
  --output-path high_quality_dataset.csv
```

### Comprehensive Kinase Dataset
```bash
capricho get \
  --target-ids CHEMBL203,CHEMBL204,CHEMBL279,CHEMBL1844 \
  --confidence-scores 7,8,9 \
  --bioactivity-type IC50,Ki,Kd \
  --calculate-pchembl \
  --aggregate-mutants \
  --max-assay-match \
  --compound-equality mixed_fp \
  --metadata-columns organism,tissue,cell_type,assay_description \
  --output-path kinase_comprehensive.csv
```

### Small Molecule Dataset via Web API
```bash
capricho get \
  --molecule-ids CHEMBL25,CHEMBL941,CHEMBL1023 \
  --chembl-backend webresource \
  --confidence-scores 8,9 \
  --output-path small_molecules.csv
```

## capricho binarize

Convert aggregated bioactivity data to binary labels (active/inactive) based on a pChEMBL threshold. This command handles censored measurements (< and >) and validates agreement between different measurement types.

```bash
capricho binarize [OPTIONS]
```

### Required Options

| Option | Description |
|---|---|
| `-i`, `--input-path` | Path to aggregated data file (CSV, TSV, or Parquet) |
| `-o`, `--output-path` | Path to save the binarized output file |

### Binarization Options

| Option | Description | Default |
|---|---|---|
| `-t`, `--threshold` | Activity threshold for binarization (pChEMBL scale) | `6.0` (1 µM) |
| `-vcol`, `--value-column` | Column to use for binarization | `pchembl_value_mean` |
| `-cid`, `--compound-id-col` | Column name for compound identifiers | `connectivity` |
| `-tid`, `--target-id-col` | Column name for target identifiers | `target_chembl_id` |
| `-rel`, `--relation-col` | Column name for standard_relation values | `standard_relation` |
| `-bcol`, `--binary-col` | Name for the output binary column | `activity_binary` |
| `-cmp-mut`, `--compare-across-mutants` | Compare measurements across mutants for conflicts | `False` |

### Understanding pChEMBL Thresholds

The pChEMBL scale is -log10(Molar), where higher values indicate higher activity:

- **pChEMBL 6.0** = 1 µM (common threshold)
- **pChEMBL 6.5** = 316 nM
- **pChEMBL 7.0** = 100 nM (stringent threshold)
- **pChEMBL 5.0** = 10 µM (permissive threshold)

### Handling Standard Relations

The binarization process handles different measurement types:

- **`=`** (discrete): Direct comparison to threshold
- **`~`** (approximate): Uses lower bound for conservative classification (±0.5 log units)
- **`<`, `<<`** (censored active): Compound is _more_ active than reported value
- **`>`, `>>`** (censored inactive): Compound is _less_ active than reported value

### Conflict Detection

The command flags measurements that disagree for the same compound-target pair:

- **Mixed discrete/censored conflicts**: When discrete measurements (`=`, `~`) disagree with censored measurements (`<`, `>`)
- **Binary label conflicts**: When measurements result in different activity classifications (active vs inactive)
- **Mutation handling**: Use `--compare-across-mutants` to control whether different mutants are compared

#### Compound Identifiers for Conflict Detection

By default, conflicts are detected using the `connectivity` column (InChI key connectivity layer), which groups compounds by their molecular graph ignoring stereochemistry. You can use a different identifier:

- **`connectivity`** (default): Groups by connectivity layer, ignoring stereochemistry
- **`smiles`**: Groups by standardized SMILES, which may be more or less permissive depending on your aggregation settings

To use SMILES for conflict detection, specify `-cid smiles`:

```bash
capricho binarize -i data.csv -o output.csv -cid smiles
```

This allows you to check for inconsistencies at different levels of molecular identity (e.g., detecting conflicts between stereoisomers when using SMILES, or treating stereoisomers as the same compound when using connectivity).

### Examples

#### Basic Binarization
```bash
# Default threshold of 6.0 (1 µM)
capricho binarize -i aggregated_data.csv -o binarized_data.csv
```

#### Custom Threshold
```bash
# Stringent threshold of 7.0 (100 nM)
capricho binarize -i egfr_data.csv -o egfr_binary.csv -t 7.0
```

#### Using Median Values
```bash
# Use median instead of mean for binarization
capricho binarize \
  -i aggregated_data.csv \
  -o binary_median.csv \
  -vcol pchembl_value_median \
  -t 6.5
```

#### Compare Across Mutants
```bash
# Flag conflicts even when measurements are on different mutants
capricho binarize \
  -i kinase_data.csv \
  -o kinase_binary.csv \
  -t 6.5 \
  --compare-across-mutants
```

### Output Format

The output file contains all original columns plus:

- **Binary activity column** (default: `activity_binary`): 0 (inactive), 1 (active), or null (missing)
- **Conflict flags**: Rows with disagreeing measurements are flagged in the `data_dropping_comment` column

Conflicting measurements are logged with detailed information about the disagreement.