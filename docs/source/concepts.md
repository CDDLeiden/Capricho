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

## Quality Control Filters

CAPRICHO provides multiple layers of quality control:

### Bioactivity Type Filtering
```bash
--bioactivity-type IC50,Ki,EC50
```
Focus on specific measurement types relevant to your analysis.

### Relation Filtering
```bash
--standard-relation "="
```
Include only exact measurements (exclude ">" or "<" relations).

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