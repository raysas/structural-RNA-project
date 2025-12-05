# WeScore: RNA Scoring Library

Goal: Creation of an objective function for the RNA folding problem

This project develops a scoring function to evaluate predicted RNA tertiary structures based on interatomic distance distributions. Supervised by Professor Guillaume Postic.

Team: Rayane Adams, Joelle Assy, Yazid Hoblos, Denys Buryi, Raul Duran De Alba

## Quick CLI Workflow (`rna-score`)

Install (editable):

```bash
pip install -r requirements.txt
pip install -e .
```

1) **Download structures**

```bash
rna-score access -n 50 -f cif -o rna_structures --workers 4
```

Could validate downloaded files with `--validate` flag.

2) **Extract distances** 

```bash
rna-score extract --folder rna_structures/mmcif --format mmcif --out-dir dist_data
```

3) **Train scoring tables**

```bash
# Histogram / simple KDE (uses scipy)
rna-score train --input-dir dist_data --out-dir training_output --method histogram

# Smoothed KDE (uses R density via rpy2 when available, falls back to scipy)
rna-score kde-train --raw-dir dist_data --output-dir kde_output --max-distance 20 --grid-step 0.1
```

4) **Score structures**

```bash
rna-score score --folder rna_structures/mmcif --tables training_output --format mmcif --output scores.csv
```

5) **Plot**

```bash
rna-score plot --input-dir training_output --output-dir plots --combined
```

Each subcommand supports `--help / -h`.
