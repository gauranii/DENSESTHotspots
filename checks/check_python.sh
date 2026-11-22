#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# File:    checks/check_python.sh
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Build every table and figure with the R scripts (scripts/) and
#          with the Python scripts (python/scripts/), time both, and compare
#          the two by code (checks/check_python.py): every cell of every csv
#          file and the marks of every figure.
# Usage:   bash checks/check_python.sh archive|rerun <work dir> [PYTHON.md]
#          [folder of the comparison tables]
#          Run from the repository root, inside `nix develop`.
#          archive: every script with its defaults: the simulation tables
#                   from the archived runs, the protein domain tables and
#                   figures from the data, the histograms of simulated data
#                   from their seed; this is what `nix flake check` runs.
#                   Outputs of the protein domain data are skipped when
#                   data/realdata.txt (or ZINULL_DATA) is not there.
#          rerun:   the 16 simulation tables by simulation, the replications
#                   RUNS (default 1:2) of every setting and true cutoff, with
#                   the seeds of the full runs.
#          The outputs go to <work dir>/<mode>/{r,python}; finished outputs
#          are kept, so an interrupted run continues. The comparison is
#          written to python_<mode>.csv in the folder given last (checks/ by
#          default), and both modes together to the Markdown file given
#          (checks/PYTHON.md by default).
#          CORES (default 4) is passed to the simulation scripts as --cores;
#          IMPLS (default "r python") chooses the implementations to run.
#          Exit status 1 if a comparison fails.
# -----------------------------------------------------------------------------
set -euo pipefail
mode=$1
work=$2
report=${3:-checks/PYTHON.md}
tables=${4:-checks}  # folder of python_<mode>.csv
cores=${CORES:-4}
impls=${IMPLS:-r python}
runs=${RUNS:-1:2}
python=${PYTHON:-python}
data=${ZINULL_DATA:-data/realdata.txt}
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# Options of one script in the given mode.
options() {
  if [ "$mode" = rerun ]; then
    echo "--rerun --runs=$runs --cores=$cores"
  fi
}

# Run a command and record its time.
timed() {
  local impl=$1 name=$2; shift 2
  local dir=$work/$mode/$impl
  mkdir -p "$dir/logs"
  [ -e "$dir/logs/$name.done" ] && return 0
  local start=$EPOCHREALTIME
  ZINULL_OUTPUT=$dir "$@" > "$dir/logs/$name.log" 2>&1 \
    || { echo "$impl $name failed, see $dir/logs/$name.log" >&2; cat "$dir/logs/$name.log" >&2; exit 1; }
  awk -v a="$start" -v b="$EPOCHREALTIME" -v c="$cores" \
    'BEGIN { printf "%.1f,%d\n", b - a, c }' > "$dir/logs/$name.done"
  echo "$impl $name: $(cat "$dir/logs/$name.done") (s, cores)"
}

for script in scripts/*.R; do
  name=$(basename "$script" .R)
  if [ "$mode" = rerun ] && ! grep -q 'simulation_table(' "$script"; then
    continue
  fi
  if [ ! -f "$data" ] && grep -q 'data/realdata.txt' "$script"; then
    echo "skipped $name ($data not found)"
    continue
  fi
  for impl in $impls; do
    if [ "$impl" = r ]; then
      timed r "$name" Rscript "scripts/$name.R" $(options)
    else
      timed python "$name" "$python" "python/scripts/$name.py" $(options)
    fi
  done
done

[ "$impls" = "r python" ] || exit 0
status=0
"$python" checks/check_python.py "$mode" "$work/$mode/r" "$work/$mode/python" \
  "$tables/python_$mode.csv" "$runs" || status=1
other=$([ "$mode" = archive ] && echo rerun || echo archive)
if [ ! -e "$tables/python_$other.csv" ] && [ -e "checks/python_$other.csv" ]; then
  cp "checks/python_$other.csv" "$tables/"
fi
"$python" checks/write_python.py "$tables/python_archive.csv" "$tables/python_rerun.csv" > "$report"
exit $status
