# Library API

* [`rna-score` Documentation](library.md)
* [User Guide](user_guide.ipynb)

```{toctree}
:hidden:
:maxdepth: 2

library.md
user_guide.ipynb
```

The `rna_score` library provides a Python API for RNA structure analysis, distance extraction, scoring, and training. This page documents the public modules, classes, and functions with their docstrings automatically extracted from the source code. In fact, this library is a wrapper around the CLI commands for programmatic access, which are based on the following scripts, whose implementation details are better described in the [CLI section](../cli.md):

```
src/
├── access_rna_structures.py
├── extract_distances.py
├── kde_training.py
├── plot_scores.py
├── score_structures.py
├── train.py
└── utils
    ├── structure_io.py
    └── validate_pdb_files.py
```

```{note}
The Kernel Density Estimation (KDE) is implemented in **R**, inside the
`kde_training.py` script, and is invoked from Python.
```

## Key Features

### Configurable Parameters

All constants in the pipeline are configurable via both CLI and API:

- **Distance cutoff** (`cutoff`): Maximum distance for contacts (default: 20.0 Å)
- **Sequence separation** (`seq_sep`): Minimum residue offset for contacts (default: 4)
- **Bin width** (`bin_width`): Histogram bin size (default: 1.0 Å)
- **Max score** (`max_score`): Score capping value (default: 10.0)
- **Pseudocount** (`pseudocount`): Smoothing parameter (default: 0.0)

### Atom Representations

The library supports multiple atomic representations:

- **Single atom** (default: `C3'`): Distance between specific atoms
- **All atoms** (`all`): Considers all heavy atoms
- **Centroid** (`centroid`): Center of mass of residues
- **Custom selection**: Specify multiple atoms (e.g., `P`, `C4'`, `O3'`)

### Scoring Formulas

Training supports multiple scoring formulas via the `scoring_formula` parameter:

| Formula | Expression | Description | Reference |
|---------|------------|-------------|----------|
| `log` | `-log(f_obs / f_ref)` | Sippl's statistical potential (default) | Standard |
| `inverse` | `f_ref / f_obs` | Inverse frequency ratio | Alternative |
| `info-gain` | `-(f_obs - f_ref) / f_ref` | Information gain (normalized) | [Postic et al., 2020](https://doi.org/10.1016/j.csbj.2020.08.013) |
| `ratio` | `f_obs / f_ref` | Direct frequency ratio | Alternative |

**Note**: The `info-gain` formula implements the **total information gain** metric from Postic et al. (2020), which is equivalent to $1 - f_{obs}/f_{ref}$ and is particularly effective for model quality assessment.

### KDE Training

The library implements Kernel Density Estimation using **R's `density()` function** via rpy2, with automatic fallback to SciPy if R is unavailable.
