# Training: from RNA structures to scoring tables

## Distance extraction

For each structure in a training set of known RNA 3D structures, we compute pairwise distances under several constraints:

- **Atom selection:** only C3′ atoms by default.
- **Intra-chain only:** distances are computed within the same chain (“intrachain”).
- **Minimum sequence separation:** only consider residue pairs separated by at least 3 positions on the sequence (e.g., $i$ with $i+4$, $i+5$, …). In the code this is implemented as a configurable threshold `--seq-sep` (default 4).
- **Cutoff:** only distances up to a maximum cutoff (default 20 Å).

`extract_distances.py` uses a KD-tree neighbor search (fast contact enumeration) and then filters pairs by same chain / same model and by sequence separation before grouping distances by base-pair type (and also pooling them into `XX`).

## Building 10 pair-specific distributions + a reference

The course baseline trains 10 distributions for the 10 unordered base pairs:

$$
\{AA, AU, AC, AG, UU, UC, UG, CC, CG, GG\}
$$

A pooled reference distribution is computed as “XX”, where residue identities are ignored.

## Histogram frequencies (course baseline)

The observed frequency in a distance bin $r$ is defined as:

$$
 f^{\mathrm{OBS}}_{i,j}(r) = \frac{N_{i,j}(r)}{N_{i,j}}
$$

The reference frequency is computed similarly over the pooled “XX” set:

$$
 f^{\mathrm{REF}}_{X,X}(r) = \frac{N_{X,X}(r)}{N_{X,X}}
$$

The spec targets 20 distance intervals between 0 and 20 Å (per pair).

## Pseudo-energy (log-odds score)

The pseudo-energy (statistical potential) for pair type $(i,j)$ at distance bin $r$ is:

$$
\bar{u}_{i,j}(r) = -\log\left(\frac{f^{\mathrm{OBS}}_{i,j}(r)}{f^{\mathrm{REF}}_{X,X}(r)}\right)
$$

To keep the potential numerically stable, the spec caps the maximum scoring value (default 10). The implementation follows this idea: `train.py` computes `-log(obs/ref)` and caps non-finite values and large values at `max_score`.

### Alternative Scoring Formulas

In addition to the standard log-odds formula (Sippl's potentials), the system supports alternative scoring functions via the `--scoring-formula` parameter:

| Formula | Expression | Description | Reference |
|---------|------------|-------------|----------|
| `log` (default) | $-kT \log\left(\frac{f^{\mathrm{OBS}}}{f^{\mathrm{REF}}}\right)$ | Sippl's statistical potential | Standard approach |
| `inverse` | $\frac{f^{\mathrm{REF}}}{f^{\mathrm{OBS}}}$ | Inverse frequency ratio (no logarithm) | Alternative metric |
| `info-gain` | $-\frac{f^{\mathrm{OBS}} - f^{\mathrm{REF}}}{f^{\mathrm{REF}}}$ | Total information gain (normalized) | Postic et al., 2020* |
| `ratio` | $\frac{f^{\mathrm{OBS}}}{f^{\mathrm{REF}}}$ | Direct frequency ratio | Enrichment metric |

**Note on `info-gain` formula**: This implements the information gain metric from [Postic et al., 2020](https://doi.org/10.1016/j.csbj.2020.08.013), which normalizes the frequency difference by the reference frequency:

$$\text{score} = -\sum_{i,j} \frac{f^{\mathrm{OBS}}_{i,j}(r) - f^{\mathrm{REF}}_{i,j}(r)}{f^{\mathrm{REF}}_{i,j}(r)}$$

These alternatives are useful for:
- Model quality assessment (distinguishing correct from incorrect structures)
- Exploring non-logarithmic scoring relationships
- Comparing different statistical frameworks
- Testing linear vs. nonlinear score transformations

Example:
```bash
# Use information gain formula (Postic et al., 2020)
rna-score train --input-dir distances --scoring-formula info-gain --output-dir tables_info_gain

# Use inverse ratio
rna-score train --input-dir distances --scoring-formula inverse --output-dir tables_inv
```

*Postic G, et al. (2020). An information gain-based approach for evaluating protein structure models. *Computational and Structural Biotechnology Journal*, 18:2228-2236. DOI: [10.1016/j.csbj.2020.08.013](https://doi.org/10.1016/j.csbj.2020.08.013)

## KDE variant (optional, smoother profiles)

Besides the histogram-based baseline, the repository also supports a KDE workflow:

- `extract_distances.py --method kde` writes raw distances per pair type.
- `kde_training.py` fits a density model on a fine grid (default step 0.1 Å) and writes frequency/score tables, with a small density floor to avoid $\log(0)$.

This typically produces smoother scoring profiles than coarse histograms.
