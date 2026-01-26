# Scoring: estimating the energy of a candidate structure

Given a candidate RNA conformation, we:

1. Compute all distances under the same constraints used in training (cutoff, intra-chain, minimum sequence separation, atom selection).
2. For each interaction distance $d$ and corresponding pair type $(i,j)$, retrieve a score $\bar{u}_{i,j}(d)$.

Because $d$ is continuous and the trained tables are discrete, the course specifies using linear interpolation to obtain a score at the exact distance.

The total score (energy estimate) for the structure is:

$$
S(\text{structure}) \,=\, \sum_{(p,q)\in\mathcal{I}} \bar{u}_{t(p),t(q)}\!\left(d(p,q)\right)
$$

where:

- $\mathcal{I}$ is the set of retained interactions (pairs of residues passing all filters),
- $t(p)$ is the nucleotide type of residue $p$,
- $d(p,q)$ is the interatomic distance (Å).

Lower scores are interpreted as more favorable conformations.

This evaluation logic is implemented in `score_structures.py`, which loads the learned `score_*.csv` tables and accumulates per-interaction scores into a total estimate.
