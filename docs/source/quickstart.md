# Quick Start Guide

This guide will help you get started with CAPRICHO in just a few minutes.

## Tab Completion

CAPRICHO supports tab completion for commands and options. To enable it, run:

```bash
capricho --install-completion
```

This will make it easier to discover available commands and options as you type.

## Basic Workflow

CAPRICHO follows a simple four-step workflow:

1. **Download** ChEMBL database (one-time setup)
2. **Explore** the data to understand what's available
3. **Get** curated bioactivity data for your analysis
4. **Binarize** (optional) convert continuous activity values to binary labels

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

## Step 4: Binarize Data (Optional)

If you need binary labels (active/inactive) for classification tasks, you can binarize the aggregated data:

### Basic Binarization

Convert continuous pChEMBL values to binary labels using the default threshold of 6.0 (1 µM):

```bash
capricho binarize -i egfr_data.csv -o egfr_binary.csv
```

### Custom Threshold

Use a more stringent threshold of 7.0 (100 nM):

```bash
capricho binarize -i egfr_data.csv -o egfr_binary.csv -t 7.0
```

### Using Median Values

Use median instead of mean for binarization:

```bash
capricho binarize -i egfr_data.csv -o egfr_binary.csv -vcol pchembl_value_median
```

The binarization process:
- Handles censored measurements (< and >) intelligently
- Flags conflicting measurements for the same compound-target pair
- Outputs a new column `activity_binary` with values 0 (inactive) or 1 (active)

## Understanding the Output

CAPRICHO generates several files:

- **Main data file** (e.g., `egfr_data.csv`): The curated bioactivity data
- **Recipe file** (e.g., `egfr_data_recipe.json`): Complete record of parameters used
- **Not aggregated file** (e.g., `egfr_data_not_aggregated.csv`): Pre-aggregation data
- **Removed subset file** (e.g., `egfr_data_removed_subset.csv`): Data filtered out during curation

### Key Output Columns

The main data file contains these important columns:

- `connectivity`: Molecular connectivity identifier
- `smiles`: Standardized SMILES representation
- `target_chembl_id`: ChEMBL ID for the target
- `pchembl_value_mean`: Mean pChEMBL value (aggregated)
- `pchembl_value_median`: Median pChEMBL value (aggregated)
- `pchembl_value_std`: Standard deviation
- `standard_relation`: Relation type (=, <, <<, >, >>, ~)
- `data_dropping_comment`: Flags for filtered data
- `data_processing_comment`: Processing notes

## Next Steps

- Read the [CLI Reference](cli-reference.md) for complete parameter documentation
- Check out [Tutorials](tutorials/index.md) for specific use cases
- Learn about [Key Concepts](concepts.md) like compound equality and backends