# DENSESTHotspots: Details

Everything beyond the overview in the [README](../README.md): the file layout, how to run and test the code, every output with its script, input and check, the data and the pinned software versions.

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
| `data/realdata.txt` | Mutation counts of the six protein domains (not in the repository, see Data) |
| `python/zinull/` | The same modules in Python (`zigp.py`, `em.py`, `cutoff.py`, `fdr.py`, `simulation.py`, `real_data.py`, `figure_data.py`, `plots.py`, `setup.py`), with `rng.py` (the random numbers of R), `rbase.py` and `rmath.py` (R's arithmetic and the probability functions it calls) |
| `python/scripts/` | One Python script per table or figure, with the names and options of `scripts/` |
| `tests/` | Unit tests of every function in `R/` (testthat): `testthat/test-<module>.R`, run by `testthat.R` |
| `python/tests/` | Unit tests of every function in `python/zinull/` (pytest): `test_<module>.py`, compared double for double with values computed by R (`make_r_values.R`, stored in `r_values.json`) |
| `checks/` | The numbers printed in the paper (`paper_values.csv`, and `extract_paper_values.py`, which read them), the comparison scripts; `check_python.sh`, which compares the outputs of the R and the Python scripts, and its result (`PYTHON.md`, with the stored comparisons `python_archive.csv` and `python_rerun.csv`); `versions.sh`, which prints the software versions |
| `flake.nix`, `flake.lock`, `scripts.nix` | Pinned environment (R with ggplot2 and testthat; Python with numpy, pandas, matplotlib, numba and pytest), one command per output in R and in Python, and the checks |
| `.github/workflows/ci.yml` | Continuous integration: `nix flake check` on every push and pull request |
| `docs/DETAILS.md` | This file |

## Running

With [Nix](https://nixos.org) (flakes enabled), from the repository root:

```sh
nix run .#table3        # one output: table1 to table3, webtableS2 to webtableS16,
                        # fig1, fig2, webfigS1 to webfigS12
nix run .#all           # every output
nix run .#py-table3     # the same output with the Python scripts (py-table1, ..., py-webfigS12, py-all)
nix flake check         # build every output, compare it with the paper, compare R with Python, run the unit tests
nix run .#test          # the unit tests of the R and the Python code
nix run .#versions      # print the software versions
nix run .#help          # list every command
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

`nix develop` opens a shell with the same R and Python, where the scripts can
be run directly: `Rscript scripts/table1_simulation_geometric.R [--rerun]
[--runs=1000] [--cores=4]` or `python python/scripts/table1_simulation_geometric.py
[--rerun]`. Each script lists its options in its header. The R and the Python
script of an output take the same options and write the same file names.

## Tests

Every function in `R/` and in `python/zinull/` has unit tests on small inputs
with a known answer (`test_every_function.py` checks that every Python
function is called by its test file).

| Suite | Files | Covers |
|---|---|---|
| R, testthat | `tests/testthat/test-<module>.R`, one per file of `R/` | ZIGP probabilities, distribution function and random numbers against their definitions and cases worked by hand; the EM steps and fits on known nulls; the cutoff criteria and the bound D; the three procedures on cases worked by hand; the simulation settings, seeds, exact storage of the replications and the first archived replication of ZIP1; the reader (with its message when the data file is missing) and the rows of Table 3 and Web Table S16; arguments and paths; the data and files of every plot |
| Python, pytest | `python/tests/test_<module>.py`, one per file of `python/zinull/` | Every function against the value R computes for the same input (`python/tests/r_values.json`, written by `make_r_values.R`), double for double; R's random numbers, arithmetic and probability functions (`test_rng.py`, `test_rbase.py`, `test_rmath.py`) |

```sh
nix run .#test                                   # both suites, in the working tree
nix develop -c Rscript tests/testthat.R          # the R suite
nix develop -c python -m pytest python/tests     # the Python suite
```

`nix flake check` runs the R suite as `checks.unit-r` and the Python suite as
`checks.unit-python`, in the sandbox, without network access.

Every function that draws random numbers takes its seed where it draws
(`generate_zigp()`, `simulate_geometric_mixture()`,
`simulate_binomial_mixture()`, and replication k of a setting with seed
889191 + k), so a result does not depend on the state of the random stream
before the call or on the number of parallel processes.

### Continuous integration

`.github/workflows/ci.yml` runs `nix flake check` on every push and pull
request (macOS on Apple silicon, Nix installed by
`cachix/install-nix-action`): the Python code matches R double for double on
arm64 macOS, where R's long double is a double (`python/zinull/rbase.py`). The
protein domain counts are not in the repository, so their outputs are skipped
there.

## Outputs

Each table and figure, the R script that builds it (the Python script has the
same name under `python/scripts/`, with `.py`), what it reads, the files it
writes (figures under `figures/`) and how `nix flake check` checks it.
Tables 1 and 2 and Web Tables S2 to S15 are built from the archived runs
(recomputed with `--rerun`); the other outputs are always recomputed.
The simulation tables give the mean and standard deviation over replications
for each null model and cutoff criterion.

| Paper | Script | Input | Output | Check |
|---|---|---|---|---|
| Table 1 | `table1_simulation_geometric.R` | `archive/simulation/` (ZIP1, ZIGP1, ZIGP2, C = 5) | `results/Table1_simulation_geometric-pi0.80_C5.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Table 2 | `table2_simulation_binomial.R` | `archive/simulation/` (ZIP2, ZIGP3, ZIGP4, C = 5) | `results/Table2_simulation_binomial-pi0.35_C5.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Table 3 | `table3_protein_rejections.R` | `data/realdata.txt` | `results/Table3_protein_rejections.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Web Tables S2 to S5 | `webtableS2_bias_geometric.R` to `webtableS5_bias_binomial.R` | `archive/simulation/` (C = 5) | `results/WebTableS<N>_bias_<f1>-pi0<pi0>_C5.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Web Tables S6, S7 | `webtableS6_simulation_geometric.R`, `webtableS7_simulation_binomial.R` | `archive/simulation/` (C = 5) | `results/WebTableS<N>_simulation_<f1>-pi0<pi0>_C5.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Web Tables S8 to S11 | `webtableS8_bias_geometric.R` to `webtableS11_bias_binomial.R` | `archive/simulation/` (C = 1) | `results/WebTableS<N>_bias_<f1>-pi0<pi0>_C1.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Web Tables S12 to S15 | `webtableS12_simulation_geometric.R` to `webtableS15_simulation_binomial.R` | `archive/simulation/` (C = 1) | `results/WebTableS<N>_simulation_<f1>-pi0<pi0>_C1.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Web Table S16 | `webtableS16_protein_estimates.R` | `data/realdata.txt` | `results/WebTableS16_protein_estimates.csv` | every cell against the paper (`checks/paper_values.csv`); R against Python, cell by cell |
| Figure 1 | `fig1_simulation_histograms.R` | none (simulated data, seed 889302) | `Fig1A_simulation_histogram_ZIP1-geometric.pdf` to `Fig1F_simulation_histogram_ZIGP4-binomial.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, mark by mark |
| Figure 2 | `fig2_protein_histograms.R` | `data/realdata.txt` | `Fig2A_protein_histogram_cd00031.pdf` to `Fig2F_protein_histogram_pfam00001.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, mark by mark |
| Web Figures S1 to S6 | `webfigS1_simulation_histograms.R` to `webfigS6_simulation_histograms.R` | none (simulated data, seed 889302) | `FigS<N><A-C>_simulation_histogram_<setting>-<f1>.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, mark by mark |
| Web Figures S7 to S12 | `webfigS7_protein_null_fit.R` to `webfigS12_protein_null_fit.R` | `data/realdata.txt` | `FigS<N><A,B>_protein_null-fit_<domain>_<C1,C2>.pdf` | marks of each PDF against the published file of the same name in `archive/figures/`; R against Python, mark by mark |

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

The comparison is written to `RESULTS.md` in the build output, one row per
table row and figure.

`checks/check_python.sh` (run by `nix flake check` as `checks.python`) builds
every output with the R and with the Python scripts and compares them: every
csv cell as text, without tolerance, and every figure by its marks (the same
bars, points and line vertices in the same order and colours, to 0.02 pt after
one straight-line map per axis). The result is `checks/PYTHON.md`: 18 of 18
csv outputs identical in every cell and 14 of 14 figures with the same marks,
and with `--rerun`, the first replications of all 24 settings identical in
every stored cell. The two agree bit for bit on Apple silicon (arm64);
x86_64 has not been compared.

## Data

- **Protein domain mutation counts** (`data/realdata.txt`): one column per
  domain, the number of somatic mutations at each position, from the
  tumors of 5848 patients in [The Cancer Genome Atlas](https://www.cancer.gov/tcga) (TCGA), mapped to
  positions within the domain: cd00031 (366 positions),
  cd00054 (90), cd00180 (871), cd00204 (409), cd00882 (762) and pfam00001
  (430). Shorter columns are padded with NA. Used by Table 3, Web Table S16,
  Figure 2 and Web Figures S7 to S12.
  The file is not included in this repository; requests go to Maricel
  Kann's laboratory (Department of Biological Sciences, University of
  Maryland, Baltimore County). With the file at
  `data/realdata.txt` of the directory the command runs from (or its path
  in `ZINULL_DATA`), these outputs build as the others do; without it,
  their scripts stop with a message and the checks skip them.
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
  numpy        2.5.2
  pandas       3.0.4
  matplotlib   3.11.1
  numba        0.67.0
  pytest       9.1.1

R        4.6.1
  boot         1.3-32
  brio         1.1.5
  callr        3.8.0
  class        7.3-23
  cli          3.6.6
  cluster      2.1.8.2
  codetools    0.2-20
  cpp11        0.5.5
  crayon       1.5.3
  desc         1.4.3
  diffobj      0.3.6
  evaluate     1.0.5
  farver       2.1.2
  foreign      0.8-91
  fs           2.1.0
  ggplot2      4.0.3
  glue         1.8.1
  gtable       0.3.6
  isoband      0.3.0
  jsonlite     2.0.0
  KernSmooth   2.23-26
  labeling     0.4.3
  lattice      0.22-9
  lifecycle    1.0.5
  magrittr     2.0.5
  MASS         7.3-65
  Matrix       1.7-5
  mgcv         1.9-4
  nlme         3.1-169
  nnet         7.3-20
  otel         0.2.0
  pkgbuild     1.4.8
  pkgload      1.5.3
  praise       1.0.0
  processx     3.9.0
  ps           1.9.3
  R6           2.6.1
  RColorBrewer 1.1-3
  rlang        1.3.0
  rpart        4.1.27
  rprojroot    2.1.1
  S7           0.2.2
  scales       1.4.0
  spatial      7.3-18
  survival     3.8-6
  testthat     3.3.2
  vctrs        0.7.3
  viridisLite  0.4.3
  waldo        0.6.2
  withr        3.0.3
```
