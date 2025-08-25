# Basic Data Retrieval

This tutorial will guide you through your first CAPRICHO query, from setup to analysis-ready data.

## Prerequisites

- CAPRICHO installed (`pip install capricho`)
- Basic familiarity with command line

## Step 1: Download ChEMBL Database

Before retrieving data, you need to download the ChEMBL database locally:

```bash
capricho download
```

This will download the latest ChEMBL version to `~/.data/chembl/`. The download is about 10GB and takes 1-2 hours, but only needs to be done once.

**Output:**
```
Downloading ChEMBL 33...
✓ Download complete: ~/.data/chembl/33/
Database ready for queries!
```

## Step 2: Explore Available Data

Let's explore what's available for a specific target. We'll use EGFR (CHEMBL203) as an example:

```bash
capricho explore --query "SELECT COUNT(*) FROM activities WHERE target_chembl_id = 'CHEMBL203'"
```

**Output:**
```
Query: SELECT COUNT(*) FROM activities WHERE target_chembl_id = 'CHEMBL203'
Results:
COUNT(*)
    45123

45,123 bioactivity measurements found for EGFR
```

## Step 3: Basic Data Retrieval

Now let's retrieve bioactivity data for EGFR:

```bash
capricho get --target-ids CHEMBL203 --output-path egfr_basic.csv
```

**What happens:**
1. CAPRICHO queries the local ChEMBL database
2. Applies default filters (confidence scores 7-9, standard bioactivity types)
3. Processes and curates the data
4. Saves results to `egfr_basic.csv`

**Output:**
```
🎯 Querying ChEMBL for target: CHEMBL203
📊 Found 45,123 raw bioactivity measurements
🔍 Applying quality filters...
   ├─ Confidence scores: 32,456 measurements (retained: 28,934)
   ├─ Bioactivity types: 28,934 measurements (retained: 24,567)
   ├─ Standard relations: 24,567 measurements (retained: 22,123)
   └─ Final dataset: 22,123 measurements
💾 Saved to: egfr_basic.csv
📋 Recipe saved: egfr_basic_recipe.json
```

## Step 4: Examine the Results

Let's look at what was generated:

```bash
ls -la egfr_basic*
```

**Output:**
```
-rw-r--r-- 1 user user 2.1M Jan 15 10:30 egfr_basic.csv
-rw-r--r-- 1 user user  1.2K Jan 15 10:30 egfr_basic_recipe.json
-rw-r--r-- 1 user user   15K Jan 15 10:30 egfr_basic.log
```

### Main Data File (egfr_basic.csv)

The main output contains these key columns:

| Column | Description |
|---|---|
| `molecule_chembl_id` | ChEMBL compound identifier |
| `target_chembl_id` | ChEMBL target identifier |
| `standard_value` | Bioactivity measurement value |
| `standard_units` | Units (usually nM) |
| `pchembl_value` | -log10 of molar potency |
| `bioactivity_type` | IC50, Ki, EC50, etc. |
| `confidence_score` | ChEMBL confidence (7-9) |

### Recipe File (egfr_basic_recipe.json)

The recipe file records exactly how the data was generated:

```json
{
  "command": "capricho get --target-ids CHEMBL203 --output-path egfr_basic.csv",
  "parameters": {
    "target_ids": ["CHEMBL203"],
    "confidence_scores": [7, 8, 9],
    "bioactivity_type": ["Potency", "Kd", "Ki", "IC50", "AC50", "EC50"],
    "standard_relation": ["="],
    "compound_equality": "connectivity"
  },
  "chembl_version": "33",
  "timestamp": "2024-01-15T10:30:25",
  "output_files": {
    "main_data": "egfr_basic.csv",
    "recipe": "egfr_basic_recipe.json",
    "log": "egfr_basic.log"
  }
}
```

## Step 5: Basic Analysis

Now you can load the data for analysis. Here's a Python example:

```python
import pandas as pd

# Load the curated data
df = pd.read_csv('egfr_basic.csv')

print(f"Dataset shape: {df.shape}")
print(f"Unique compounds: {df['molecule_chembl_id'].nunique()}")
print(f"Bioactivity types: {df['bioactivity_type'].value_counts()}")
print(f"Confidence distribution: {df['confidence_score'].value_counts()}")
```

**Output:**
```
Dataset shape: (22123, 15)
Unique compounds: 8956
Bioactivity types: 
IC50        12456
Ki           5678
EC50         2345
Kd           1644
Potency       789
AC50          211

Confidence distribution:
9    15678
8     4567
7     1878
```

## Next Steps

Congratulations! You've successfully retrieved and examined your first CAPRICHO dataset. Here's what to explore next:

1. **Quality Control**: Learn about [filtering options](quality-control.md) to refine your data
2. **Understanding Output**: Deep dive into [interpreting results](understanding-output.md)
3. **Target Analysis**: Follow the [target analysis tutorial](target-analysis.md) for biological insights

## Common Issues

**"No data found" error**: Check that the target ID is correct and exists in ChEMBL
**Large file warnings**: For targets with >50k measurements, consider additional filtering
**Memory issues**: Use `--skip-not-aggregated` to save space for large datasets

## Summary

In this tutorial, you:
- Downloaded the ChEMBL database
- Explored available data
- Retrieved your first curated dataset
- Examined the output structure
- Started basic analysis

The key insight is that CAPRICHO provides not just the data, but complete transparency about how it was processed, ensuring your analysis is reproducible and trustworthy.