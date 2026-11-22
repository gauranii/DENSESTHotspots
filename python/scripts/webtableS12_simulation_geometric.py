# -----------------------------------------------------------------------------
# File:    python/scripts/webtableS12_simulation_geometric.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Web Table S12.
# Purpose: Number of rejections, FDR and TPR of the two-stage, one-stage and
#          Storey procedures (mean and standard deviation over replications)
#          for the settings ZIP1, ZIGP1, ZIGP2: geometric alternative, pi0 = 0.80,
#          true cutoff C = 1 (contaminated null), four null models, two cutoff criteria.
# Inputs:  archive/simulation/<setting>_C1.csv.gz (all replications of the
#          three settings); none with --rerun (simulated data).
# Outputs: results/WebTableS12_simulation_geometric-pi0.80_C1.csv; with --rerun also
#          results/simulation/<setting>_C1/runNNNN.csv.
# Usage:   python python/scripts/webtableS12_simulation_geometric.py [--rerun] [--runs=1000] [--cores=4]
#          --rerun: run the replications here instead of reading the archive.
#          --runs: N (runs 1 to N), a:b or a comma list; finished runs are kept.
#          Options and outputs as scripts/webtableS12_simulation_geometric.R, which this file follows;
#          --help prints this header.
# Source:  scripts/webtableS12_simulation_geometric.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PAPER, PATHS, parse_args  # noqa: E402
from zinull.simulation import simulation_table  # noqa: E402

opts = parse_args({"runs": str(PAPER["replications"]), "cores": 4})

simulation_table("WebTableS12", "simulation", f1="geometric", pi0=0.80, C=1, opts=opts,
                 settings=PAPER, paths=PATHS)
