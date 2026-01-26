# Outputs and how they are used

- **Training outputs** (`training_output/` or a user-specified directory):
  - `freq_XX.csv`, `freq_<PAIR>.csv`: estimated distance frequencies.
  - `score_<PAIR>.csv`: the potential $\bar{u}_{i,j}(r)$ as a function of distance.
  - `metadata.json`: parameters used during training.
- **Plotting outputs** (`plot_scores.py`):
  - per-pair scoring profiles (static PNG and/or interactive HTML).
- **Structure scoring outputs** (`score_structures.py`):
  - total score per structure, plus optional per-interaction details.

Workflow summary: known structures → distance statistics → statistical potential → scoring of candidate folds.
