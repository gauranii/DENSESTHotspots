# R and Python

Written by `checks/check_python.sh`: the outputs of the R scripts (`scripts/`)
and of the Python scripts (`python/scripts/`) compared by code
(`checks/check_python.py`). Every csv file under `results/` is compared cell
by cell, without tolerance: the same header, the same rows, every cell the
same text, and the same bytes. Every figure is compared by the marks read from
the two PDF files (`checks/figure_compare.py`): the same number of bars, or of
points and line vertices in each series, in the same order and with the same
colours; bar ends, bar heights, point centres and line vertices agree to
0.02 pt of the R page after one straight-line map per axis (one scale for the
bar heights), since R and matplotlib lay their pages out differently and R
prints its coordinates with two decimals. "Values" counts the compared cells
and mark coordinates; "Largest difference" is in the units of the cells, or in
pt for a figure. Times are wall-clock seconds on an 8 core Apple M2 laptop
shared with other jobs; the Python times include loading its compiled (numba)
functions.

The tests in `python/tests/` (also run by `nix flake check`) compare every
function of `python/zinull/` with values computed by R
(`python/tests/make_r_values.R`), double for double.

## With the defaults of the scripts

Both implementations run every script with its defaults, as `nix flake check`
runs them (`bash checks/check_python.sh archive <work dir>`): the 16
simulation tables from the archived runs (1000 replications per setting), the
histograms of simulated data from their seed, and Table 3, Web Table S16,
Figure 2 and Web Figures S7 to S12 from the protein domain counts
(`data/realdata.txt`; these are skipped when the file is not there).

18 of 18 outputs in csv files identical in every cell; 14 of 14 figures with the same marks (4 cores).

| Output | From | Files | Values | Not identical | Largest difference | Verdict | R (s) | Python (s) |
|---|---|---|---|---|---|---|---|---|
| Table 1 | archived runs | 1 | 624 | 0 | 0 | identical | 0.7 | 0.8 |
| Table 2 | archived runs | 1 | 624 | 0 | 0 | identical | 0.8 | 0.7 |
| Table 3 | computed from data | 1 | 252 | 0 | 0 | identical | 7.1 | 3.6 |
| Figure 1 | seeded simulated data | 6 | 1068 | 0 | 0.00726 pt | same marks | 0.3 | 2.7 |
| Figure 2 | computed from data | 6 | 1347 | 0 | 0.00747 pt | same marks | 0.2 | 1.0 |
| Web Table S2 | archived runs | 1 | 432 | 0 | 0 | identical | 0.7 | 0.7 |
| Web Table S3 | archived runs | 1 | 432 | 0 | 0 | identical | 0.8 | 0.7 |
| Web Table S4 | archived runs | 1 | 432 | 0 | 0 | identical | 0.9 | 0.7 |
| Web Table S5 | archived runs | 1 | 432 | 0 | 0 | identical | 0.9 | 0.7 |
| Web Table S6 | archived runs | 1 | 624 | 0 | 0 | identical | 0.9 | 0.7 |
| Web Table S7 | archived runs | 1 | 624 | 0 | 0 | identical | 0.9 | 0.7 |
| Web Table S8 | archived runs | 1 | 432 | 0 | 0 | identical | 0.8 | 0.7 |
| Web Table S9 | archived runs | 1 | 432 | 0 | 0 | identical | 0.8 | 0.7 |
| Web Table S10 | archived runs | 1 | 432 | 0 | 0 | identical | 0.9 | 0.7 |
| Web Table S11 | archived runs | 1 | 432 | 0 | 0 | identical | 1.0 | 0.7 |
| Web Table S12 | archived runs | 1 | 624 | 0 | 0 | identical | 0.8 | 0.7 |
| Web Table S13 | archived runs | 1 | 624 | 0 | 0 | identical | 0.8 | 0.7 |
| Web Table S14 | archived runs | 1 | 624 | 0 | 0 | identical | 0.9 | 0.7 |
| Web Table S15 | archived runs | 1 | 624 | 0 | 0 | identical | 1.0 | 0.7 |
| Web Table S16 | computed from data | 1 | 432 | 0 | 0 | identical | 7.1 | 1.8 |
| Web Figure S1 | seeded simulated data | 3 | 594 | 0 | 0.00552 pt | same marks | 0.2 | 1.0 |
| Web Figure S2 | seeded simulated data | 3 | 435 | 0 | 0.00679 pt | same marks | 0.2 | 1.0 |
| Web Figure S3 | seeded simulated data | 3 | 561 | 0 | 0.00652 pt | same marks | 0.2 | 1.0 |
| Web Figure S4 | seeded simulated data | 3 | 462 | 0 | 0.00679 pt | same marks | 0.2 | 1.0 |
| Web Figure S5 | seeded simulated data | 3 | 549 | 0 | 0.00722 pt | same marks | 0.2 | 1.0 |
| Web Figure S6 | seeded simulated data | 3 | 435 | 0 | 0.00504 pt | same marks | 0.2 | 1.0 |
| Web Figure S7 | computed from data | 2 | 124 | 0 | 0.00498 pt | same marks | 1.5 | 1.0 |
| Web Figure S8 | computed from data | 2 | 144 | 0 | 0.00469 pt | same marks | 1.6 | 1.0 |
| Web Figure S9 | computed from data | 2 | 132 | 0 | 0.00531 pt | same marks | 1.5 | 1.0 |
| Web Figure S10 | computed from data | 2 | 116 | 0 | 0.00535 pt | same marks | 1.3 | 1.0 |
| Web Figure S11 | computed from data | 2 | 100 | 0 | 0.00589 pt | same marks | 1.2 | 1.0 |
| Web Figure S12 | computed from data | 2 | 372 | 0 | 0.00521 pt | same marks | 4.5 | 1.9 |

## Recomputed by simulation

Both implementations run the 16 simulation tables with `--rerun` for the
first replications of all 24 combinations of a setting and a true cutoff,
with the seeds of the full runs (`bash checks/check_python.sh rerun <work
dir>`). Compared are the stored file of every replication (rejections, FDR
and TPR of the three procedures, estimates and their bias, with 17
significant digits) and the two tables of each group of settings. The times
are those of the two table scripts of a group together; the second one reads
the replications the first one stored.

8 of 8 outputs in csv files identical in every cell (4 cores).

| Output | From | Files | Values | Not identical | Largest difference | Verdict | R (s) | Python (s) |
|---|---|---|---|---|---|---|---|---|
| Table 1, Web Table S2 | geometric, pi0 = 0.80, C = 5: replications 1:2 of ZIP1, ZIGP1, ZIGP2 | 8 | 2160 | 0 | 0 | identical | 3.3 | 1.9 |
| Table 2, Web Table S3 | binomial, pi0 = 0.35, C = 5: replications 1:2 of ZIP2, ZIGP3, ZIGP4 | 8 | 2160 | 0 | 0 | identical | 3.5 | 1.4 |
| Web Table S4, Web Table S6 | geometric, pi0 = 0.85, C = 5: replications 1:2 of ZIP3, ZIGP5, ZIGP6 | 8 | 2160 | 0 | 0 | identical | 3.0 | 1.4 |
| Web Table S5, Web Table S7 | binomial, pi0 = 0.80, C = 5: replications 1:2 of ZIP4, ZIGP7, ZIGP8 | 8 | 2160 | 0 | 0 | identical | 2.5 | 1.3 |
| Web Table S8, Web Table S12 | geometric, pi0 = 0.80, C = 1: replications 1:2 of ZIP1, ZIGP1, ZIGP2 | 8 | 2160 | 0 | 0 | identical | 4.0 | 1.5 |
| Web Table S9, Web Table S13 | binomial, pi0 = 0.35, C = 1: replications 1:2 of ZIP2, ZIGP3, ZIGP4 | 8 | 2160 | 0 | 0 | identical | 2.8 | 1.3 |
| Web Table S10, Web Table S14 | geometric, pi0 = 0.85, C = 1: replications 1:2 of ZIP3, ZIGP5, ZIGP6 | 8 | 2160 | 0 | 0 | identical | 3.1 | 1.4 |
| Web Table S11, Web Table S15 | binomial, pi0 = 0.80, C = 1: replications 1:2 of ZIP4, ZIGP7, ZIGP8 | 8 | 2160 | 0 | 0 | identical | 2.5 | 1.3 |

