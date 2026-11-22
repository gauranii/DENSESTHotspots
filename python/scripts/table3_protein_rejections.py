# -----------------------------------------------------------------------------
# File:    python/scripts/table3_protein_rejections.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471. Table 3.
# Purpose: Numbers of positions rejected at FDR level 0.10 in six protein domains by
#          the one-stage, two-stage and Storey procedures, for four null models and
#          two cutoff criteria.
# Inputs:  data/realdata.txt.
# Outputs: results/Table3_protein_rejections.csv.
# Usage:   python python/scripts/table3_protein_rejections.py
#          Options and outputs as scripts/table3_protein_rejections.R, which this file follows;
#          --help prints this header.
# Source:  scripts/table3_protein_rejections.R.
# -----------------------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from zinull.setup import PAPER, PATHS, output_file, show_help  # noqa: E402

show_help()

from zinull.rbase import write_csv  # noqa: E402
from zinull.real_data import protein_rejections  # noqa: E402

table = protein_rejections(PATHS["data"], PAPER["fdr_level_real"])
write_csv(table, output_file(PATHS["results"], "Table3_protein_rejections.csv"))
print(table.to_string(index=False))
