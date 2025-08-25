# Quick Start Guide

This guide will help you get started with CAPRICHO in just a few minutes.

## Basic Workflow

CAPRICHO follows a simple three-step workflow:

1. **Download** ChEMBL database (one-time setup)
2. **Explore** the data to understand what's available
3. **Get** curated bioactivity data for your analysis

## Step 1: Download ChEMBL Database

First, download the ChEMBL database to your local machine:

```bash
capricho download
```

This downloads the latest ChEMBL version to `~/.data/chembl/` by default. You can specify a different version or location:

```bash
capricho download --version 33 --prefix /path/to/custom/location
```

## Step 2: Explore the Database

Before fetching data, you might want to explore what's available:

```bash
# List all tables in the database
capricho explore --list-tables

# Explore a specific table
capricho explore --table activities

# Search for tables containing a specific column
capricho explore --search-column pchembl

# Run a custom SQL query
capricho explore --query "SELECT COUNT(*) FROM activities WHERE confidence_score = 9"
```

## Step 3: Get Bioactivity Data

Now you're ready to fetch and curate bioactivity data. Here are some common use cases:

### Basic Target Query

Get all bioactivity data for EGFR (CHEMBL203):

```bash
capricho get --target-ids CHEMBL203 --output-path egfr_data.csv
```

### Multiple Targets with High Confidence

Get data for multiple targets, but only high-confidence measurements:

```bash
capricho get --target-ids CHEMBL203,CHEMBL204,CHEMBL279 --confidence-scores 8,9 --output-path high_confidence.csv
```

### Specific Bioactivity Types

Focus on specific bioactivity measurements:

```bash
capricho get --target-ids CHEMBL203 --bioactivity-type IC50,Ki --output-path ic50_ki_data.csv
```

### With Custom Aggregation

Aggregate mutant data and include additional metadata:

```bash
capricho get --target-ids CHEMBL203 --aggregate-mutants --metadata-columns organism,tissue --output-path aggregated_data.csv
```

## Understanding the Output

CAPRICHO generates several files:

- **Main data file** (e.g., `egfr_data.csv`): The curated bioactivity data
- **Recipe file** (e.g., `egfr_data_recipe.json`): Complete record of parameters used
- **Log file**: Detailed processing information

### Key Output Columns

The main data file contains these important columns:

- `molecule_chembl_id`: ChEMBL ID for the compound
- `target_chembl_id`: ChEMBL ID for the target
- `standard_value`: Bioactivity measurement
- `standard_units`: Units of measurement
- `pchembl_value`: Calculated or reported pChEMBL value
- `confidence_score`: ChEMBL confidence score (0-9)
- Various metadata columns

## Next Steps

- Read the [CLI Reference](cli-reference.md) for complete parameter documentation
- Check out [Tutorials](tutorials/index.md) for specific use cases
- Learn about [Key Concepts](concepts.md) like compound equality and backends