# EmpiricalNullZIGP: empirical null estimation with zero-inflated discrete mixtures

R code for

> Gauran, I. I. M., Park, J., Lim, J., Park, D., Zylstra, J., Peterson, T.,
> Kann, M., and Spouge, J. L. (2018). Empirical null estimation using
> zero-inflated discrete mixture distributions and its application to protein
> domain data. *Biometrics* 74:458-471. doi:10.1111/biom.12779.

The method finds the unusually large counts in a vector of counts, for
example the positions of a protein domain that carry more somatic mutations
in tumors than the rest. The counts at or below an unknown cutoff C follow the null;
the null is a zero-inflated generalized Poisson (ZIGP), zero-inflated Poisson
(ZIP), generalized Poisson (GP) or Poisson (P) distribution, fitted by an EM
algorithm that treats the counts above C as missing. C is chosen by one of two
criteria (C1, C2), and the positions are then tested by three procedures: the
one-stage procedure (local FDR), the two-stage procedure (local FDR below a
bound D, rejection at or above D) and Storey's FDR procedure. Each table and
figure of the paper and of its web supplement has one script in `scripts/`,
and `nix flake check` builds all of them and compares them with the paper.
Figure 3 of the paper shows protein structures and has no code.

```r
source("R/setup.R")
x <- read_domains(PATHS$data)$cd00054        # mutation counts of 90 positions
fits <- select_null_fits(x, cap = max(x))    # four null models, two cutoff criteria
fits
two_stage_count(x, fits[1, ], level = 0.10)  # positions rejected: ZIGP null, criterion C1
```

## Layout

| Path | Content |
|---|---|
| `R/setup.R` | Loader sourced by every script: settings of the paper, paths, argument parser |
| `R/zigp.R` | ZIGP random numbers, probability functions, the null mass below a cutoff and the bound D |
| `R/em.R` | EM algorithms of the four null models, fitted below every candidate cutoff |
| `R/cutoff.R` | The cutoff criteria C1 and C2, the estimate of pi0 and the bound D |
| `R/fdr.R` | One-stage, two-stage and Storey procedures, for simulated data (with labels) and real data |
| `R/simulation.R` | The 12 simulation settings, one replication, replications run in parallel and stored, and the tables of means and standard deviations |
| `R/real_data.R` | Reader of the protein domain counts, Table 3 and Web Table S16 |
| `R/figure_data.R` | Data of the figures: the simulated data sets of the histograms and the observed and fitted null probabilities |
| `R/plots.R` | Histograms and null-fit plots |
| `scripts/` | One script per table or figure |
| `archive/simulation/` | All replications of the simulation study, one file per setting and true cutoff, which the table scripts read |
| `archive/figures/` | The published figure files that come from code, under the names the scripts give their figures |
| `archive/collect_simulation.R` | Combines the replication files of a rerun into `archive/simulation/` |
| `data/realdata.txt` | Mutation counts of the six protein domains |
| `checks/` | The numbers printed in the paper (`paper_values.csv`, and `extract_paper_values.py`, which read them), the comparison scripts, and `versions.sh`, which prints the software versions |
| `flake.nix`, `flake.lock`, `scripts.nix` | Pinned environment (R with ggplot2, Python for the checks), one command per output, and the checks |

## Running

With [Nix](https://nixos.org) (flakes enabled), from the repository root:

```sh
nix run .#table3        # one output: table1 to table3, webtableS2 to webtableS16,
                        # fig1, fig2, webfigS1 to webfigS12
nix run .#all           # every output
nix flake check         # build every output and compare it with the paper
nix run .#versions      # print the software versions
```

Tables are written to `results/` (csv) and figures to `figures/` (pdf), each
under the name of its table or figure in the paper (see Outputs).

Table 3, Web Table S16 and all figures are always computed, from the data in
`data/realdata.txt` or from simulated data. Tables 1 and 2 and Web Tables S2
to S15 summarize a simulation study of 1000 replications for each of 12 null
settings and two true cutoffs (C = 5, and C = 1 for the contaminated null);
by default their scripts read the stored replications in
`archive/simulation/`. To run the replications here instead, add `--rerun`:

```sh
nix run .#table1 -- --rerun --runs=1:100 --cores=6
```

Replication k of a setting uses the seed 889191 + k, so the results do not
depend on the number of cores or on the order of the runs. With `--rerun`
each replication is stored as `results/simulation/<setting>_C<C>/runNNNN.csv`;
finished replications are kept, so a study can be run in pieces, and the two
tables of a group of settings (for example Table 1 and Web Table S2) share
their replications. `Rscript archive/collect_simulation.R` then combines them
into `archive/simulation/`.

`nix develop` opens a shell with the same R, where the scripts can be run
directly: `Rscript scripts/table1_simulation_geometric.R [--rerun]
[--runs=1000] [--cores=4]`. Each script lists its options in its header.

## Outputs

Each table and figure, the script that builds it and the files it writes.
The simulation tables give the mean and standard deviation over replications
for each null model and cutoff criterion.

| Paper | Script | Input | Output |
|---|---|---|---|
| Table 1 | `table1_simulation_geometric.R` | `archive/simulation/` (ZIP1, ZIGP1, ZIGP2, C = 5) | `results/Table1_simulation_geometric-pi0.80_C5.csv` |
| Table 2 | `table2_simulation_binomial.R` | `archive/simulation/` (ZIP2, ZIGP3, ZIGP4, C = 5) | `results/Table2_simulation_binomial-pi0.35_C5.csv` |
| Table 3 | `table3_protein_rejections.R` | `data/realdata.txt` | `results/Table3_protein_rejections.csv` |
| Web Tables S2 to S5 | `webtableS2_bias_geometric.R` to `webtableS5_bias_binomial.R` | `archive/simulation/` (C = 5) | `results/WebTableS<N>_bias_<f1>-pi0<pi0>_C5.csv` |
| Web Tables S6, S7 | `webtableS6_simulation_geometric.R`, `webtableS7_simulation_binomial.R` | `archive/simulation/` (C = 5) | `results/WebTableS<N>_simulation_<f1>-pi0<pi0>_C5.csv` |
| Web Tables S8 to S11 | `webtableS8_bias_geometric.R` to `webtableS11_bias_binomial.R` | `archive/simulation/` (C = 1) | `results/WebTableS<N>_bias_<f1>-pi0<pi0>_C1.csv` |
| Web Tables S12 to S15 | `webtableS12_simulation_geometric.R` to `webtableS15_simulation_binomial.R` | `archive/simulation/` (C = 1) | `results/WebTableS<N>_simulation_<f1>-pi0<pi0>_C1.csv` |
| Web Table S16 | `webtableS16_protein_estimates.R` | `data/realdata.txt` | `results/WebTableS16_protein_estimates.csv` |
| Figure 1 | `fig1_simulation_histograms.R` | none (simulated data, seed 889302) | `Fig1A_simulation_histogram_ZIP1-geometric.pdf` to `Fig1F_simulation_histogram_ZIGP4-binomial.pdf` |
| Figure 2 | `fig2_protein_histograms.R` | `data/realdata.txt` | `Fig2A_protein_histogram_cd00031.pdf` to `Fig2F_protein_histogram_pfam00001.pdf` |
| Web Figures S1 to S6 | `webfigS1_simulation_histograms.R` to `webfigS6_simulation_histograms.R` | none (simulated data, seed 889302) | `WebFigS<N><A-C>_simulation_histogram_<setting>-<f1>.pdf` |
| Web Figures S7 to S12 | `webfigS7_protein_null_fit.R` to `webfigS12_protein_null_fit.R` | `data/realdata.txt` | `WebFigS<N><A,B>_protein_null-fit_<domain>_<C1,C2>.pdf` |

Commands and run times:

| Paper | Command | Default | With `--rerun` |
|---|---|---|---|
| Tables 1, 2, Web Tables S2 to S15 | `nix run .#table1` and so on | about 1 s each | 18 to 30 min per group of three settings, 6 cores |
| Table 3, Web Table S16 | `nix run .#table3`, `nix run .#webtableS16` | about 20 s each (always computed from the data) | same |
| Figures and web figures | `nix run .#fig1` and so on | a few seconds each (always computed) | same |
| All 24 simulation settings | the 8 simulation table scripts with `--rerun` | | about 3 h, 6 cores |

Times are from an 8 core Apple laptop shared with other work.

## Checks

`checks/run_checks.sh` (run by `nix flake check`) builds every output with
the defaults of its script and compares it with the paper:

- every cell of Tables 1 to 3 and Web Tables S2 to S16 against the number
  printed in the paper (`checks/paper_values.csv`), after rounding the
  computed value to the decimals of the printed cell (`check_tables.py`);
- every figure against the published file of the same name in
  `archive/figures/`, by the marks in the two PDF files
  (`check_figures.py`): the bars of each histogram (same number of bars,
  heights proportional), and the points of each series of a null-fit plot
  (same number of points, same positions after one affine map per axis).

The comparison is written to `checks/RESULTS.md`, one row per table row
and figure.

## Data

- **Protein domain mutation counts** (`data/realdata.txt`): one column per
  domain, the number of somatic mutations at each position, from the
  tumors of 5848 patients in The Cancer Genome Atlas (TCGA), mapped to
  positions within the domain: cd00031 (366 positions),
  cd00054 (90), cd00180 (871), cd00204 (409), cd00882 (762) and pfam00001
  (430). Shorter columns are padded with NA. Used by Table 3, Web Table S16,
  Figure 2 and Web Figures S7 to S12.
- **Stored replications** (`archive/simulation/<setting>_C<C>.csv.gz`): for
  every replication, cutoff criterion and null model, the number of
  rejections R, the FDR and the TPR of the three procedures, the estimates
  and their bias. The FDR of a replication without rejections is NaN, and
  the tables average the FDR over the replications with rejections. Numbers are written with 17 significant digits and read
  back exactly. ZIP1 to ZIGP8 are the 12 settings of the paper; the tables
  of the contaminated null (C = 1) label them MZIP1 to MZIGP8.
- **Published figures** (`archive/figures/`): the figure files of the
  article and the web supplement that come from code.

## Software versions

`flake.lock` pins the nixpkgs revision, and with it the version of R, Python
and every package (there is no `renv.lock`). `nix run .#versions` (or
`bash checks/versions.sh` inside `nix develop`) prints them; `nix flake check`
compares this section with its output.

```text
nixpkgs  NixOS/nixpkgs f3aa7d3d72c2bc66e0f30f4262d7ec2e33c9eab2

Python   3.14.7

R        4.6.1
  boot         1.3-32
  class        7.3-23
  cli          3.6.6
  cluster      2.1.8.2
  codetools    0.2-20
  cpp11        0.5.5
  farver       2.1.2
  foreign      0.8-91
  ggplot2      4.0.3
  glue         1.8.1
  gtable       0.3.6
  isoband      0.3.0
  KernSmooth   2.23-26
  labeling     0.4.3
  lattice      0.22-9
  lifecycle    1.0.5
  MASS         7.3-65
  Matrix       1.7-5
  mgcv         1.9-4
  nlme         3.1-169
  nnet         7.3-20
  R6           2.6.1
  RColorBrewer 1.1-3
  rlang        1.3.0
  rpart        4.1.27
  S7           0.2.2
  scales       1.4.0
  spatial      7.3-18
  survival     3.8-6
  vctrs        0.7.3
  viridisLite  0.4.3
  withr        3.0.3
```
