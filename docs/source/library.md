# Library API

The `rna_score` library provides a Python API for RNA structure analysis, distance extraction, scoring, and training. This page documents the public modules, classes, and functions with their docstrings automatically extracted from the source code.

## Quick Start

```python
from rna_score import RNAScorer

# Initialize the scorer
scorer = RNAScorer()

# Extract distances from RNA structures
scorer.extract_distances(folder="rna_structures/pdb", out_dir="distances")

# Train scoring tables
scorer.train_scoring(output_dir="training_output")

# Plot scoring profiles
scorer.plot_scores(output_dir="plots", combined=True)

# Score a structure
scorer.score_structure(pdb_path="test_structure.pdb", output_csv="scores.csv")
```

## Module Functions

```{eval-rst}
.. autofunction:: rna_score.lib.download_rna_structures
```

## RNAScorer Class

```{eval-rst}
.. autoclass:: rna_score.lib.RNAScorer
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__, __call__
```
