# DENSESTHotspots: Discrete Empirical Null for Sparse Estimation and Simultaneous Testing

Code for

> Gauran, I. I. M., Park, J., Lim, J., Park, D., Zylstra, J., Peterson, T.,
> Kann, M., and Spouge, J. L. (2018). Empirical Null Estimation Using
> Zero-Inflated Discrete Mixture Distributions and its Application to Protein
> Domain Data. *Biometrics* 74:458-471.
> [doi.org/10.1111/biom.12779](https://doi.org/10.1111/biom.12779)

## The idea

Which positions of a protein domain carry more somatic mutations in tumors
than chance would give? Most positions carry few or none, so the counts have
far more zeros than a Poisson null allows, and the null that is chosen
decides how many positions are found. The paper takes the counts at or below an unknown cutoff C as the
null and fits a zero-inflated generalized Poisson (ZIGP), zero-inflated
Poisson (ZIP), generalized Poisson (GP) or Poisson (P) distribution to them by
an EM algorithm that treats the counts above C as missing. C is chosen by one
of two criteria (C1, C2), and the positions are then tested by a one-stage
(local FDR), a two-stage and Storey's FDR procedure. The paper studies this
in simulations and on six protein domains from The Cancer Genome Atlas.

```r
source("R/setup.R")
x <- read_domains(PATHS$data)$cd00054        # mutation counts of 90 positions
fits <- select_null_fits(x, cap = max(x))    # four null models, two cutoff criteria
two_stage_count(x, fits[1, ], level = 0.10)  # positions rejected: ZIGP null, criterion C1
```

## What this repository does

It rebuilds every table and figure of the paper and its web supplement that
comes from code (Tables 1 to 3, Web Tables S2 to S16, Figures 1 and 2, Web
Figures S1 to S12; Figure 3 shows protein structures and has no code) and
checks each one against the published paper. The simulation tables read the
archived replications; everything else is recomputed.

The code exists twice, in R (`R/`, `scripts/`) and in Python (`python/`).
The Python version is an independent second implementation, a check of the
paper's numbers: it gives the same numbers as R, cell by cell for the tables
and mark by mark for the figures, so the results do not depend on one
language or one set of libraries. A Nix flake pins R, Python and every
package, so the checks run the same way on any machine.

## Running

With [Nix](https://nixos.org) (flakes enabled), from the repository root:

```sh
nix run .#table3        # one output: table1 to table3, webtableS2 to webtableS16, fig1, fig2, webfigS1 to webfigS12
nix run .#all           # every output, in R
nix run .#py-table3     # the same output from the Python code (py-all for every one)
nix flake check         # build everything, compare with the paper and R with Python, run the tests
nix run .#help          # every command
```

Tables go to `results/` (csv) and figures to `figures/` (pdf), named after
the paper. `checks/PYTHON.md` records the comparison of R with Python. Every
function in R and Python has unit tests, and GitHub runs `nix flake check` on
every push.

## Data

- **Protein domain mutation counts** (`data/realdata.txt`): somatic mutations
  at each position of six domains, from the tumors of 5848 patients in
  [The Cancer Genome Atlas](https://www.cancer.gov/tcga). The file is not in
  the repository; requests go to Maricel Kann's laboratory (Department of
  Biological Sciences, University of Maryland, Baltimore County). Without it,
  Table 3, Web Table S16, Figure 2 and Web Figures S7 to S12 are skipped.
- **Simulation study**: every replication is archived in
  `archive/simulation/`; `--rerun` recomputes them.

## More detail

[`docs/DETAILS.md`](docs/DETAILS.md) has the file layout, every output with
its script, input and check, run times, the tests, the data and the pinned
software versions.

## Citation

Please cite the paper above when you use this code.
