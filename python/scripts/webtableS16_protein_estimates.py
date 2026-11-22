# -----------------------------------------------------------------------------
# File:    python/scripts/webtableS16_protein_estimates.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Web Table S16.
# Purpose: Estimates of eta, lambda, theta, pi0, the cutoff C and the bound D in six
#          protein domains, for four null models and two cutoff criteria.
# Inputs:  data/realdata.txt.
# Outputs: results/WebTableS16_protein_estimates.csv.
# Usage:   python python/scripts/webtableS16_protein_estimates.py
#          Options and outputs as scripts/webtableS16_protein_estimates.R, which this file follows;
#          --help prints this header.
# Source:  scripts/webtableS16_protein_estimates.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PATHS, output_file, show_help  # noqa: E402

show_help()

from zinull.rbase import write_csv  # noqa: E402
from zinull.real_data import protein_estimates  # noqa: E402

table = protein_estimates(PATHS["data"])
write_csv(table, output_file(PATHS["results"], "WebTableS16_protein_estimates.csv"))
print(table.to_string(index=False))
