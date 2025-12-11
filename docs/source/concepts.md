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

### SMILES-Based

The `smiles` method uses standardized SMILES strings directly for exact matching:

```bash
capricho get --target-ids CHEMBL203 --compound-equality smiles
```

**Advantages:**
- Simple and transparent matching logic
- No additional computation (InChI or fingerprints)
- Useful when you trust the standardized SMILES

**Note:** This method relies on the standardization performed by the ChEMBL structure pipeline. Beware that different tautomers may not match even if they represent the same compound. Connectivities can be more robust, but may merge stereoisomers.

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

## Non-pChEMBL Aggregation

By default, CAPRICHO aggregates bioactivity data using the `pchembl_value` column, which represents -log10(molar) potency values. However, many ChEMBL assays (especially ADMET assays) report measurements that aren't suitable for pChEMBL conversion, such as:

- **Permeability** (e.g., Caco-2 apparent permeability in cm/s)
- **Percent inhibition** (e.g., % inhibition at a fixed concentration)
- **Clearance** (e.g., mL/min/kg)
- **Half-life** (e.g., hours)
- **Solubility** (e.g., µg/mL)

For these measurements, use `--aggregate-on standard_value`:

```bash
capricho get \
  --assay-ids CHEMBL1112933,CHEMBL3529279 \
  --assay-types A \
  --aggregate-on standard_value \
  --output-path permeability_data.csv
```

### Key Differences from pChEMBL Aggregation

| Aspect | pchembl_value (default) | standard_value |
|--------|------------------------|----------------|
| **Mean type** | Geometric mean | Arithmetic mean |
| **Units** | Always molar (-log10) | Original units preserved |
| **Use case** | Potency measurements (IC50, Ki, etc.) | ADMET, physicochemical properties |

### Important: Preventing Unit Mixing

When aggregating on `standard_value`, ensure you don't inadvertently combine measurements with different units. Use `--id-columns` to group by unit type:

```bash
capricho get \
  --assay-ids CHEMBL1112933,CHEMBL3529279 \
  --aggregate-on standard_value \
  --id-columns standard_units \
  --output-path permeability_data.csv
```

This creates separate aggregations for measurements with different `standard_units` values.

## Unit Conversion

ChEMBL contains bioactivity data reported in many different units, even for the same measurement type. For example, permeability might be reported as `cm/s`, `nm/s`, or `10^-6 cm/s`. This heterogeneity makes cross-study aggregation challenging.

The `--convert-units` flag enables automatic unit conversion before aggregation:

```bash
capricho get \
  --assay-ids CHEMBL1112933,CHEMBL3529279 \
  --aggregate-on standard_value \
  --convert-units \
  --output-path permeability_data.csv
```

### Supported Unit Families

| Family | Target Unit | Source Units |
|--------|------------|--------------|
| **Permeability** | `10^-6 cm/s` | `cm/s`, `nm/s`, `ucm/s`, `10'-6 cm/s`, etc. |
| **Molar concentration** | `nM` | `uM`, `µM`, `mM`, `pM`, `M` |
| **Mass concentration** | `ug/mL` | `ng/ml`, `mg/ml`, `mg/L`, `pg/ml` |
| **Dose** | `mg/kg` | `ug/kg`, `ug.kg-1`, `mg.kg-1` |
| **Time** | `hr` | `min`, `s`, `ms`, `day` |

### Transparency

All unit conversions are logged and tracked in the `data_processing_comment` column. This ensures you can always trace which measurements were converted and by what factor.

### Example: Caco-2 Permeability Dataset

```bash
capricho get \
  --assay-ids CHEMBL1112933,CHEMBL3529279,CHEMBL3529278 \
  --assay-types A \
  --confidence-scores 0,1,2,3,4,5,6,7,8,9 \
  --aggregate-on standard_value \
  --convert-units \
  --id-columns standard_units,assay_cell_type \
  --drop-unassigned-chiral \
  --output-path caco2_permeability.csv
```

This command:
1. Fetches permeability data from specific Caco-2 assays
2. Converts all permeability units to `10^-6 cm/s`
3. Aggregates using arithmetic mean on `standard_value`
4. Groups by cell type to maintain biological context

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

### pchembl_relation: Inverted Relations for -log Scale

When working with pChEMBL values ($-log_{10}(Molar)$), the direction of comparison operators is inverted compared to the original concentration values. CAPRICHO automatically creates a `pchembl_relation` column during binarization to make this relationship explicit:

**Relation Inversion Logic:**
- `standard_relation` `<` (low concentration, active) → `pchembl_relation` `>` (high pChEMBL, active)
- `standard_relation` `>` (high concentration, inactive) → `pchembl_relation` `<` (low pChEMBL, inactive)
- `standard_relation` `=` → `pchembl_relation` `=` (unchanged)
- `standard_relation` `~` → `pchembl_relation` `~` (unchanged)

**Example Interpretation:**
For a measurement with `IC50` = 1 µM (pChEMBL = 6.0) and `standard_relation` = `<`:
- Original: IC50 < 1 µM (active at concentrations below 1 µM)
- pChEMBL: pchembl_value > 6.0 (higher pChEMBL = more active)
- With threshold = 6.0: classified as **active (1)**

This inverted relation column is automatically added when you run the `binarize` command and helps interpret how measurements relate to activity thresholds on the -log scale.

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