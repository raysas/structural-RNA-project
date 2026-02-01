# RNA folding as an optimization problem

For a given RNA chain (a sequence of ribonucleotides), the RNA folding problem consists of finding the native 3D conformation among an astronomically large number of possible structures. In thermodynamic terms, the native fold is typically modeled as the conformation with the lowest Gibbs free energy.

Because computing an exact physical free energy for an arbitrary RNA conformation is difficult, practical methods often rely on an objective function that approximates (or correlates with) the Gibbs free energy. In this project, the objective function is a knowledge-based statistical potential learned from experimentally solved RNA structures.

If a geometric pattern appears frequently in experimentally determined RNA structures, it is assumed to be energetically favorable. Conversely, rare patterns are treated as less favorable. Patterns are modeled through interatomic distances between nucleotides (here: using only the C3′ atom by default, or a configurable alternative).

## Configurable Parameters

The pipeline provides flexibility through several configurable parameters:

- **Atom representation**: C3′ (default), all atoms, centroid, or specific atom selections
- **Distance cutoff**: Maximum contact distance (default: 20 Å)
- **Sequence separation**: Minimum residue offset for contacts (default: 4)
- **Bin width**: Histogram discretization (default: 1.0 Å)
- **Scoring formula**: Log-odds (default), inverse ratio, info-gain, or direct ratio

## Statistical Potential

For each nucleotide-type pair $ (i, j) \in \{A, U, C, G\}^2 $, we estimate:

- an observed distance distribution $ f^{\mathrm{OBS}}_{i,j}(r) $ from real structures,
- a reference distribution $ f^{\mathrm{REF}}_{X,X}(r) $ (where residue identity is ignored).

A pseudo-energy is then derived from the log-ratio of these distributions, so that distances enriched in real data obtain lower scores.
