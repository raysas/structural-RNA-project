# Command-Line Interface (CLI)

The `rna-score` CLI provides a comprehensive toolkit for downloading, processing, and scoring RNA 3D structures using statistical models of interatomic distances.

## Installation

Install the package using git, as mentioned in [main page's installation instructions](./index.rst):


```bash
pip install git+https://github.com/raysas/structural-RNA-project.git
```

or clone it and install:
```bash
git clone https://github.com/raysas/structural-RNA-project.git
cd structural-RNA-project
pip install -r requirements.txt
pip install -e .
```

## Quick Start

```bash
# -- download RNA structures
rna-score access -n 50 --rna-only -o rna_structures

# -- extract distances
rna-score extract --folder rna_structures/pdb --out-dir distances

# -- train scoring tables
rna-score train --input-dir distances --output-dir training_output

# -- score a structure
rna-score score --pdb structure.pdb --tables training_output --output scores.csv

# -- visualize scoring profiles
rna-score plot --input-dir training_output --output-dir plots --combined
```

## Available Commands

The `rna-score` CLI includes the following subcommands:

- `access` - Download RNA structures from RCSB PDB
- `extract` - Extract interatomic distances from structures
- `train` - Train scoring tables from distance distributions
- `score` - Score RNA structures using trained tables
- `plot` - Visualize scoring profiles
- `workflow` - Run the complete pipeline (extract → train → score)

Use `rna-score <command> --help` to see detailed options for each command.

---

## `rna-score access`

Download RNA structures from the RCSB PDB database.

### Syntax

```bash
rna-score access [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-n, --max-structures` | int | 10 | Maximum number of structures to download |
| `--download-all` | flag | - | Download all available RNA structures |
| `--rna-only` | flag | - | Restrict to RNA-only structures (no protein/DNA) |
| `-f, --formats` | list | pdb, cif | File formats to download (pdb, cif) |
| `-o, --output-dir` | str | rna_structures | Output directory for structures |
| `--workers` | int | 5 | Number of parallel download workers |
| `--show-info` | flag | - | Display metadata for structures before downloading |
| `--list-only` | flag | - | Only retrieve PDB IDs without downloading files |
| `--validate` | flag | - | Validate downloaded files and filter out invalid ones |

### Examples

Download 50 RNA-only structures in mmCIF format:

```bash
rna-score access -n 50 --rna-only -f cif -o data/rna_structures
```

Download all RNA structures with validation:

```bash
rna-score access --download-all --rna-only --validate -o rna_structures
```

List available RNA structures without downloading:

```bash
rna-score access -n 100 --list-only --show-info
```

---

## `rna-score extract`

Extract interatomic distances from RNA structures and generate distance distributions.

### Syntax

```bash
rna-score extract [INPUT] [OPTIONS]
```

### Input Options (choose one)

| Option | Description |
|--------|-------------|
| `--pdb PDB_ID` | Download and process a single PDB ID |
| `--pdb-list FILE` | Process structures from a list file |
| `--folder DIR` | Process all structures in a directory |

**List file format** (for `--pdb-list`):
```
<PDB_ID> [CHAIN_ID1 CHAIN_ID2 ...]
```

Example:
```
1EHZ A
1Y26 B C
2OEU
```

### Processing Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--format {pdb,mmcif}` | str | pdb | Structure file format (auto-detected for local files) |
| `--atom-mode` | list | C3' | Atom selection: `C3'`, `centroid`, `all`, or specific atoms (e.g., `P C4`) |
| `--dist-mode {intra,inter}` | str | intra | Distance mode: `intra` (within chain) or `inter` (between chains) |
| `--cutoff` | float | 20.0 | Maximum distance (Å) to consider as contact |
| `--seq-sep` | int | 4 | Minimum sequence separation for intra-chain contacts (residues) |
| `--bin-width` | float | 1.0 | Bin width (Å) for histogram discretization |
| `--method {histogram,kde}` | str | histogram | Extraction method: binned counts or raw distances for KDE |
| `--cores` | int | all | Number of CPU cores for parallel processing |
| `--all-models` | flag | - | Process all NMR models (default: first model only) |
| `--save-detailed` | flag | - | Save detailed CSV with all distance measurements |
| `--out-dir` | str | dist_data | Output directory for distance distributions |

### Features Summary

| Feature | Description | Option | Details |
|---------|-------------|--------|---------|
| **Input Source** | Structure selection | `--pdb`, `--pdb-list`, `--folder` | Choose one: PDB ID, list file, or directory |
| **Format** | File format | `--format {pdb,mmcif}` | Default: pdb (auto-detected for local files) |
| **Atomic Selection** | Structure representation | `--atom-mode` | C3' (default), centroid, all, or specific atoms |
| **Interaction Mode** | Distance measurement | `--dist-mode {intra,inter}` | intra: within chain; inter: between chains |
| **Sequence Separation** | Minimum residue offset | `--seq-sep` | Default: 4 residues (ignored in inter mode) |
| **Distance Cutoff** | Maximum contact distance | `--cutoff` | Default: 20.0 Å |
| **Output Type** | Distribution format | `--method {histogram,kde}` | histogram: binned counts; kde: raw distances |
| **Parallelization** | CPU cores | `--cores` | Default: all available |
| **NMR Models** | Model processing | `--all-models` | Default: first model only |
| **Detailed Log** | Full distance info | `--save-detailed` | Saves PDB, chains, residues, atoms, distances, etc. |

### Examples

Extract distances from a folder of mmCIF files:

```bash
rna-score extract --folder rna_structures/mmcif --format mmcif --out-dir dist_data
```

Extract using multiple atom types:

```bash
rna-score extract --folder structures/ --atom-mode P C4 --method kde
```

Extract inter-chain distances with detailed output:

```bash
rna-score extract --pdb-list structures.txt --dist-mode inter --save-detailed
```

Process a single PDB structure:

```bash
rna-score extract --pdb 1EHZ --chains A --out-dir distances
```

---

## `rna-score train`

Train scoring tables from extracted distance distributions using statistical potentials or KDE.

### Syntax

```bash
rna-score train [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input-dir` | str | dist_data | Directory with distance distributions |
| `--output-dir` | str | training_output | Output directory for scoring tables |
| `--method {histogram,kde}` | str | histogram | Training method |
| `--scoring-formula {log,inverse,difference,ratio}` | str | log | Scoring formula (see below) |
| `--max-score` | float | 10.0 | Maximum score cap to prevent outliers |
| `--pseudocount` | float | 0.0 | Pseudocount for smoothing (prevents division by zero) |
| `--cutoff` | float | 20.0 | Maximum distance (Å) - must match extraction |
| `--bin-width` | float | 1.0 | Bin width (Å) - must match extraction |

### Scoring Formulas

The `--scoring-formula` parameter controls how scores are computed from observed and reference frequency distributions:

| Formula | Expression | Description | Use Case | Reference |
|---------|------------|-------------|----------|----------|
| **log** (default) | `-log(f_obs / f_ref)` | Sippl's statistical potential | Default choice for knowledge-based scoring | Standard |
| **inverse** | `f_ref / f_obs` | Inverse frequency ratio | Linear relationship without logarithm | Alternative |
| **info-gain** | `-(f_obs - f_ref) / f_ref` | Information gain (normalized) | Model quality assessment | Postic+ 2020 |
| **ratio** | `f_obs / f_ref` | Direct frequency ratio | Enrichment/depletion scores | Alternative |

**Notes:**
- The `log` formula is the standard Sippl's potential used in most knowledge-based potentials in structural biology
- The `info-gain` formula implements the **total information gain** approach from [Postic et al., 2020](https://doi.org/10.1016/j.csbj.2020.08.013), particularly effective for distinguishing correct from incorrect structural models
- Alternative formulas are useful for exploring different scoring schemes or when non-logarithmic relationships are expected
- The formula choice is saved in `metadata.json` for reproducibility

### Examples

Train scoring tables using default histogram method:

```bash
rna-score train --input-dir dist_data --output-dir training_output
```

Train with smoothing using pseudocount:

```bash
rna-score train --input-dir dist_data --output-dir tables --pseudocount 1.0 --max-score 15.0
```

Train using KDE method:

```bash
rna-score train --input-dir kde_distances --output-dir kde_tables --method kde
```

Train with alternative scoring formula:

```bash
# Standard log formula
rna-score train --input-dir dist_data --output-dir tables_log --scoring-formula log

# Inverse ratio formula
rna-score train --input-dir dist_data --output-dir tables_inverse --scoring-formula inverse

# Information gain formula
rna-score train --input-dir dist_data --output-dir tables_infogain --scoring-formula info-gain
```

### Output Files

The training process generates:
- `score_<PAIR>.csv` - Scoring tables for each base pair type (AA, AC, AG, AU, CC, CG, CU, GG, GU, UU)
- `metadata.json` - Training parameters and configuration

---

## `rna-score score`

Score RNA structures using trained scoring tables.

### Syntax

```bash
rna-score score [INPUT] [OPTIONS]
```

### Input Options (choose one)

| Option | Description |
|--------|-------------|
| `--pdb FILE` | Score a single PDB/mmCIF file |
| `--folder DIR` | Score all structures in a directory |
| `--list FILE` | Score structures from a list file |

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--tables` | str | training_output | Directory containing scoring tables |
| `--format {pdb,mmcif}` | str | pdb | Input file format |
| `--atom-mode` | list | C3' | Atom selection (must match training) |
| `--cutoff` | float | 20.0 | Distance cutoff (must match training) |
| `--seq-sep` | int | 4 | Sequence separation (must match training) |
| `--detailed` | flag | - | Print detailed per-interaction scores |
| `--output` | str | - | Output CSV file for scores |

### Examples

Score a single structure:

```bash
rna-score score --pdb test_structure.pdb --tables training_output --output scores.csv
```

Score multiple structures from a folder:

```bash
rna-score score --folder structures/ --tables training_output --format mmcif --output batch_scores.csv
```

Score with detailed output:

```bash
rna-score score --pdb structure.pdb --tables tables/ --detailed --output detailed_scores.csv
```

### Output Format

The output CSV contains:
- Structure ID
- Total score
- Per-pair-type scores (optional)
- Number of contacts (optional)

---

## `rna-score plot`

Generate visualization plots of scoring profiles.

### Syntax

```bash
rna-score plot [OPTIONS]
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--input-dir` | str | training_output | Directory with scoring tables |
| `--output-dir` | str | plots | Output directory for plots |
| `--combined` | flag | - | Generate combined multi-panel plot |

### Examples

Plot individual scoring profiles:

```bash
rna-score plot --input-dir training_output --output-dir plots
```

Generate combined visualization:

```bash
rna-score plot --input-dir training_output --output-dir plots --combined
```

### Output Files

The plotting generates:
- `png/` - Individual PNG plots for each pair type
- `html/` - Interactive Plotly HTML plots
- Combined plots (if `--combined` is used)

---

## `rna-score workflow`

Run the complete pipeline: extract distances → train scoring tables → score structures.

### Syntax

```bash
rna-score workflow [OPTIONS]
```

### Options

| Option | Type | Description |
|--------|------|-------------|
| `--train-folder` | str | Directory with training structures |
| `--train-list` | str | List file with training structures |
| `--score-folder` | str | Directory with structures to score |
| `--score-list` | str | List file with structures to score |
| `--output-dir` | str | Base output directory (default: workflow_output) |
| `--format` | str | File format (default: pdb) |
| `--method` | str | Training method (default: histogram) |
| `--scoring-formula` | str | Scoring formula (default: log) |
| `--atom-mode` | list | Atom selection (default: C3') |
| `--cutoff` | float | Distance cutoff (default: 20.0) |
| `--seq-sep` | int | Sequence separation (default: 4) |
| `--cores` | int | CPU cores for parallel processing |

### Examples

Full workflow with training and scoring:

```bash
# Create scoring list
cat <<EOF > tests/scoring_list.txt
1EHZ A
1Y26 B C
EOF

rna-score workflow \
  --train-folder rna_structures/mmcif \
  --score-list tests/scoring_list.txt \
  --output-dir workflow_output \
  --format mmcif \
  --method histogram \
  --scoring-formula log
```

Workflow with alternative scoring formula:

```bash
rna-score workflow \
  --train-folder structures/ \
  --score-list to_score.txt \
  --scoring-formula info-gain \
  --output-dir results_infogain
```

Workflow with custom parameters:

```bash
rna-score workflow \
  --train-folder structures/ \
  --score-list to_score.txt \
  --atom-mode P C4 \
  --cutoff 25.0 \
  --cores 8 \
  --scoring-formula info-gain \
  --output-dir results
```

### Output Structure

The workflow creates:
```
workflow_output/
├── extracted/          # Distance distributions
├── training_output/    # Scoring tables
└── scores.csv         # Final scores
```

---

## Getting Help

For detailed help on any command:

```bash
rna-score --help
rna-score <command> --help
```

For version information:

```bash
rna-score --version
```

## Common Workflows

### Basic Workflow

```bash
# 1. Download structures
rna-score access -n 100 --rna-only -o rna_structures

# 2. Extract distances
rna-score extract --folder rna_structures/pdb --out-dir distances

# 3. Train scoring tables
rna-score train --input-dir distances --output-dir tables

# 4. Score structures
rna-score score --folder test_structures/ --tables tables --output scores.csv

# 5. Visualize
rna-score plot --input-dir tables --output-dir plots --combined
```

### KDE-Based Workflow

```bash
# Extract raw distances for KDE
rna-score extract --folder structures/ --method kde --out-dir kde_data

# Train using KDE
rna-score train --input-dir kde_data --method kde --output-dir kde_tables

# Score with KDE tables
rna-score score --pdb structure.pdb --tables kde_tables --output scores.csv
```

### Inter-Chain Analysis

```bash
# Extract inter-chain distances
rna-score extract --folder structures/ --dist-mode inter --out-dir inter_distances

# Train and score
rna-score train --input-dir inter_distances --output-dir inter_tables
rna-score score --pdb complex.pdb --tables inter_tables --output inter_scores.csv
```

## Tips

- **Consistency**: Always use the same `--cutoff`, `--seq-sep`, and `--atom-mode` values across extraction, training, and scoring steps.
- **Parallelization**: Use `--cores` to speed up processing of large datasets.
- **Validation**: Use `--validate` with `access` command to ensure structure quality.
- **Detailed Logs**: Use `--save-detailed` during extraction for debugging or analysis of individual contacts.
- **Format Detection**: For local files, format is auto-detected; specify `--format` explicitly for downloaded PDB IDs.