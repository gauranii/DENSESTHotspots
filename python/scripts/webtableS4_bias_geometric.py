# -----------------------------------------------------------------------------
# File:    python/scripts/webtableS4_bias_geometric.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Web Table S4.
# Purpose: Bias of the estimates of eta, lambda, theta, pi0 and the cutoff C
#          (mean and standard deviation over replications)
#          for the settings ZIP3, ZIGP5, ZIGP6: geometric alternative, pi0 = 0.85,
#          true cutoff C = 5, four null models, two cutoff criteria.
# Inputs:  archive/simulation/<setting>_C5.csv.gz (all replications of the
#          three settings); none with --rerun (simulated data).
# Outputs: results/WebTableS4_bias_geometric-pi0.85_C5.csv; with --rerun also
#          results/simulation/<setting>_C5/runNNNN.csv.
# Usage:   python python/scripts/webtableS4_bias_geometric.py [--rerun] [--runs=1000] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          Options and outputs as in the R script of the same name;
#          --help prints this header.
# Source:  scripts/webtableS4_bias_geometric.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PAPER, PATHS, parse_args  # noqa: E402
from zinull.simulation import simulation_table  # noqa: E402

opts = parse_args({"runs": str(PAPER["replications"]), "cores": 4})

simulation_table("WebTableS4", "bias", f1="geometric", pi0=0.85, C=5, opts=opts,
                 settings=PAPER, paths=PATHS)
