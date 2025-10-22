# Key Concepts

Understanding these core concepts will help you use CAPRICHO effectively and make informed decisions about data curation.

## Compound Equality

One of the most important decisions in bioactivity data analysis is determining when two compound entries represent the same molecule.

### Connectivity-Based (Default)

The `connectivity` method considers compounds the same if they have identical molecular connectivity, ignoring stereochemistry:

```bash
capricho get --target-ids CHEMBL203 --compound-equality connectivity
```

**Advantages:**
- Aggregates stereoisomers together
- Useful when stereochemistry data is inconsistent
- Larger datasets with more statistical power

**Use When:**
- Stereochemistry information is unreliable
- You want to combine data from different stereoisomers
- Building large-scale datasets

### Fingerprint-Based

The `mixed_fp` method uses molecular fingerprints (ECFP4 + RDKit) to determine similarity:

```bash
capricho get --target-ids CHEMBL203 --compound-equality mixed_fp
```

**Advantages:**
- More chemically precise
- Maintains stereochemical distinctions
- Better for SAR analysis

**Use When:**
- Stereochemistry is important for your analysis
- Building focused datasets for SAR studies
- Working with well-curated data

## ChEMBL Backends

CAPRICHO supports two different ways to access ChEMBL data, each with distinct advantages.

### Local Database Backend (Default)

The `downloader` backend uses a local SQLite database:

```bash
capricho download  # One-time setup
capricho get --target-ids CHEMBL203 --chembl-backend downloader
```

**Advantages:**
- Much faster for large queries
- Works offline after initial download
- Consistent performance
- Full SQL query capabilities

**Requirements:**
- ~25GB disk space for full ChEMBL database
- Initial download time

**Best For:**
- Large-scale data mining
- Repeated queries
- Complex filtering requirements
- Offline analysis

### Web API Backend

The `webresource` backend queries the live ChEMBL web API:

```bash
capricho get --target-ids CHEMBL203 --chembl-backend webresource
```

**Advantages:**
- No local storage required
- Always uses latest data
- No setup time
- Good for small queries

**Limitations:**
- Slower for large queries
- Requires internet connection
- Subject to API rate limits

**Best For:**
- Small, targeted queries
- One-off analyses
- When disk space is limited

## Data Aggregation

CAPRICHO provides several options for handling duplicate measurements and aggregating data.

### Target Mutations

By default, CAPRICHO treats different target mutations as separate entities. Use `--aggregate-mutants` to combine them:

```bash
capricho get --target-ids CHEMBL203 --aggregate-mutants
```

This is useful when you want to study the target in general rather than specific mutations.

### Assay Matching

The `--max-assay-match` option enforces strict metadata matching between assays:

```bash
capricho get --target-ids CHEMBL203 --max-assay-match
```

This ensures that only assays with identical metadata (organism, tissue, etc.) are aggregated together.

### Metadata Columns

Include additional metadata in your analysis:

```bash
capricho get --target-ids CHEMBL203 --metadata-columns organism,tissue,cell_type
```

These columns are preserved during aggregation and can help you understand data heterogeneity.

## Confidence Scoring

ChEMBL assigns confidence scores (0-9) based on target specificity:

### High Confidence (8-9)
- **9**: Direct single protein target
- **8**: Direct protein complex/family target

Best for focused target analysis and SAR studies.

### Medium Confidence (6-7)
- **7**: Direct target with homologues
- **6**: Molecular target assigned

Good balance of data quantity and quality.

### Lower Confidence (0-5)
- Include broader, less specific data
- Useful for large-scale analyses
- May include off-target effects

## Standard Relations and Censored Data

ChEMBL bioactivity measurements include a `standard_relation` field that indicates the relationship between the measured value and the reported concentration.

### Relation Types

- **`=`**: The measured value equals the reported concentration
- **`<`**: The compound is active at concentrations _below_ the reported value (_active_ at the reported concentration)
- **`<<`**: Stronger indication that the compound is active at concentrations _well below_ the reported value (_active_ at the reported concentration)
- **`>`**: The compound is active at concentrations _above_ the reported value (_inactive_ at the reported concentration)
- **`>>`**: Stronger indication that the compound is active at concentrations _well above_ the reported value (_inactive_ at the reported concentration)
- **`~`**: Approximate measurement (CAPRICHO handles these as ±0.5 log units)

### Working with Censored Data

**Important**: ChEMBL only pre-calculates pchembl_value for exact measurements (`standard_relation='='`). To include censored data (`<`, `>`), you *must* use the `--calculate-pchembl` flag:

```bash
# Default: Only includes exact measurements (=)
capricho get --target-ids CHEMBL203 --output-path egfr_exact.csv

# Include censored measurements: MUST use --calculate-pchembl
capricho get --target-ids CHEMBL203 \
  --standard-relation "=,<,>" \
  --calculate-pchembl \
  --output-path egfr_all.csv
```

Without `--calculate-pchembl`, you'll get an error if you request censored data, but the data will still be fetched:
```
ERROR: pchembl_values are only calculated for standard_relation='='.
If you want to use censored data, please set calculate_pchembl to True.
```

### Aggregation with Censored Data

When aggregating data with censored measurements, CAPRICHO only combines measurements that have:
1. Identical `standard_relation` values
2. In case of identical (`=`) standard_relation, statistics will be calculated for compound-target pairs with multiple exact measurements.
3. Censored measurements (`<`, `>`, `<<`, `>>`) are only combined with exact measurement matches (e.g.: < 6.0 will not be combined with < 5.0).

This conservative approach prevents mixing incompatible measurement types (e.g., averaging an exact value with a lower bound).

### Activity Data Analysis

**For binary classification** (active/inactive), use the `binarize` command which properly handles censored measurements. Following the example above, we have:

```bash
capricho binarize -i egfr_all.csv -o egfr_binary.csv -t 6.0
```
See the CLI reference for detailed binarization options.

## Quality Control Filters

CAPRICHO provides multiple layers of quality control:

### Bioactivity Type Filtering
```bash
--bioactivity-type IC50,Ki,EC50
```
Focus on specific measurement types relevant to your analysis.

### Relation Filtering
```bash
# Only exact measurements (default behavior, pchembl pre-calculated by ChEMBL)
--standard-relation "="

# Include censored data (requires --calculate-pchembl)
--standard-relation "=,<,>" --calculate-pchembl
```
Choose measurement precision level based on your analysis needs.

### Date Requirements
```bash
--require-doc-date
```
Ensure all data has associated publication dates.

### Assay Size Constraints
```bash
--min-assay-size 10 --max-assay-size 1000
```
Filter assays by number of tested compounds.

## Reproducibility

CAPRICHO ensures full reproducibility through several mechanisms:

### Recipe Files

Every run generates a JSON recipe file containing:
- All command-line parameters
- ChEMBL version used
- Timestamp and environment info
- Data processing steps

```json
{
  "command": "capricho get --target-ids CHEMBL203",
  "parameters": {
    "target_ids": ["CHEMBL203"],
    "confidence_scores": [7, 8, 9],
    "bioactivity_type": ["IC50", "Ki"]
  },
  "chembl_version": "33",
  "timestamp": "2024-01-15T10:30:00"
}
```

### Version Control

Specify exact ChEMBL versions:
```bash
capricho get --target-ids CHEMBL203 --chembl-version 33
```

### Transparent Processing

All filtering steps are logged and flagged data is preserved for inspection.

## Output Structure

Understanding CAPRICHO's output helps you make the most of your curated data:

### Main Data File
- Curated bioactivity measurements
- Standardized column names
- Quality flags and metadata

### Pre-aggregation Data
- Raw data before aggregation steps
- Useful for understanding curation impact
- Can be skipped with `--skip-not-aggregated`

### Log Files
- Detailed processing information
- Quality control statistics
- Warning and error messages

This multi-file approach ensures transparency while providing clean, analysis-ready data.