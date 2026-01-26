# Command Reference

## Installation

Install (editable mode for development):

```bash
    pip install -r requirements.txt
    pip install -e .
```

## 1. Download RNA Structures

```bash
    rna-score access -n 50 --rna-only -f cif -o data/rna_structures --workers 4
```

*Add --validate to filter out invalid downloaded files.*

## 2. Extract Interatomic Distances

```bash
    rna-score extract --folder rna_structures/mmcif --format mmcif --out-dir dist_data
```

## 3. Train Scoring Tables

```bash
    rna-score train --input-dir dist_data --output-dir training_output --method histogram
```
## 4. Score Structures
```bash
    rna-score score --folder rna_structures/mmcif --tables training_output --format mmcif --output scores.csv
```
## 5. Plot Results
```bash
    rna-score plot --input-dir training_output --output-dir plots --combined
```
## 6. Full Workflow (All Steps in One Command)

First, create a scoring list with PDB IDs and chains:
```
    cat <<EOF > tests/scoring_list.txt
    1EHZ A
    1Y26 B C
    EOF
```

Then run the full workflow:
```bash
    rna-score workflow \
        --train-folder data/rna_structures/mmcif \
        --score-list tests/scoring_list.txt \
        --output-dir tests/workflow_output \
        --format mmcif \
        --method histogram
```
This command will run extraction, training, scoring, and plotting in one go.

## 7. Help for Each Subcommand

Every subcommand supports --help or -h for detailed instructions. Examples:
```bash
    rna-score access --help
    rna-score train --help
    rna-score plot --help
    rna-score extract --help
    rna-score workflow --help
    rna-score score --help
```