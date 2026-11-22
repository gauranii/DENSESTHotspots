#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/versions.sh
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Print the versions of the software the flake provides: the nixpkgs
#          revision of flake.lock, Python and the packages of pyEnv in
#          flake.nix, and R and its packages.
# Usage:   bash checks/versions.sh            print the versions
#          bash checks/versions.sh --check    compare them with the section
#                                             "Software versions" of
#                                             docs/DETAILS.md
#          Run from the repository root, inside `nix develop`.
# Outputs: prints the versions; with --check exit status 1 if docs/DETAILS.md
#          differs.
# -----------------------------------------------------------------------------
set -euo pipefail
python=${PYTHON:-python3}

versions() {
  "$python" - <<'PY'
import json
from importlib.metadata import version
import platform

locked = json.load(open("flake.lock"))["nodes"]["nixpkgs"]["locked"]
print(f"nixpkgs  {locked['owner']}/{locked['repo']} {locked['rev']}")
print()
print(f"Python   {platform.python_version()}")
for name in ["numpy", "pandas", "matplotlib", "numba", "pytest"]:  # pyEnv of flake.nix
    print(f"  {name:12} {version(name)}")
PY
  echo
  Rscript -e '
    cat(sprintf("R        %s.%s\n", R.version$major, R.version$minor))
    installed <- installed.packages()
    installed <- installed[is.na(installed[, "Priority"]) | installed[, "Priority"] != "base", ]
    installed <- installed[order(tolower(rownames(installed))), ]
    cat(sprintf("  %-12s %s\n", rownames(installed), installed[, "Version"]), sep = "")'
}

if [ "${1:-}" = --check ]; then
  diff <(versions) <(awk '/^## Software versions/ { section = 1 }
                          section && /^```/ { block++; next }
                          section && block == 1 { print }' docs/DETAILS.md)
else
  versions
fi
