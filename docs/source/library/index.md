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