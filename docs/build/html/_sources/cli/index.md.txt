# Command-Line Interface (CLI)

The `rna-score` CLI provides a comprehensive toolkit for downloading, processing, and scoring RNA 3D structures using statistical models of interatomic distances.

* [CLI Usage](usage.md)
* [Implementation Details](implementation.md)

```{toctree}
:hidden:
:maxdepth: 2

usage.md
implementation.md
```

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
