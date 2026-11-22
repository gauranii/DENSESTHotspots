# -----------------------------------------------------------------------------
# File:    python/tests/test_simulation.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of simulation.py: the simulated data, one replication, the
#          stored files and the tables, against R: the same doubles and the
#          same files, line for line.
# Usage:   python -m pytest python/tests/test_simulation.py
# -----------------------------------------------------------------------------
import gzip
import os

import numpy as np
import pandas as pd
import pytest

from conftest import numbers, same
from inputs import BINOMIAL, SETTING, SMALL
from zinull import rbase, simulation
from zinull.setup import REPO_ROOT


def test_settings():
    settings = simulation.SIMULATION_SETTINGS
    assert len(settings) == 12 and settings["setting"].is_unique
    assert list(settings.groupby(["f1", "pi0"]).size()) == [3, 3, 3, 3]
    assert simulation.F1_BINOMIAL == {"size": 250, "prob": 0.20}
    assert simulation.F1_GEOMETRIC == {"prob": 0.08}
    assert len(simulation.TABLE_COLUMNS["simulation"]) == 9
    assert simulation.TABLE_COLUMNS["bias"][0] == "bias_eta"


def test_simulate_labelled_counts(r):
    data = simulation.simulate_labelled_counts(60, 5, SETTING, 7)
    same(np.concatenate((data["value"], data["label"])), r("simulation_counts"))
    assert (data["value"][data["label"] == 1] > 5).all()  # non-null counts are above C
    data = simulation.simulate_labelled_counts(60, 1, BINOMIAL, 8)
    same(np.concatenate((data["value"], data["label"])), r("simulation_counts_binomial"))


def test_simulation_replication_and_its_file(r, tmp_path):
    replication = simulation.simulation_replication(SETTING, 5, 3, SMALL)
    assert len(replication) == 8 and list(replication.columns[:3]) == ["run", "cutoff", "model"]
    same(numbers(replication), r("simulation_replication"))
    simulation.write_exact_csv(replication, tmp_path / "run0003.csv")
    assert (tmp_path / "run0003.csv").read_text().splitlines() \
        == r.lines("simulation_replication_file")


def test_format_exact_and_parse_exact(r):
    assert simulation.format_exact(0.1) == "0.10000000000000001"
    assert simulation.format_exact(12.0) == "12"
    assert simulation.format_exact(rbase.NA) == "NA"
    assert simulation.format_exact(float("nan")) == "NaN"
    assert simulation.format_exact(float("-inf")) == "-Inf"
    parsed = simulation.parse_exact(["0.10000000000000001", "0.30000000000000004", "NA", "12"])
    same(parsed, r("simulation_parse_exact"))
    assert parsed[0] == 0.1 and parsed[1] == 0.1 + 0.2 and rbase.is_na(parsed[2])


def test_write_exact_csv_and_read_exact_csv(tmp_path):
    table = pd.DataFrame({"run": [2, 1], "cutoff": ["C1", "C2"], "model": ["P", "GP"],
                          "eta": [1 / 3, rbase.NA], "FDR": [np.nan, 0.5]})
    simulation.write_exact_csv(table, tmp_path / "t.csv")
    assert (tmp_path / "t.csv").read_text() == ("run,cutoff,model,eta,FDR\n"
                                                "2,C1,P,0.33333333333333331,NaN\n"
                                                "1,C2,GP,NA,0.5\n")
    read = simulation.read_exact_csv(tmp_path / "t.csv")
    assert list(read["run"]) == [1, 2] and list(read["model"]) == ["GP", "P"]  # sorted by run
    assert read["eta"][1] == 1 / 3 and rbase.is_na(read["eta"][0])


@pytest.fixture(scope="module")
def four_runs(tmp_path_factory):
    """Four replications of ZIGP3 with true cutoff 1, stored and read back
    in another order, as in make_r_values.R."""
    folder = tmp_path_factory.mktemp("runs")
    for run in range(1, 5):
        replication = simulation.simulation_replication(BINOMIAL, 1, run, SMALL)
        simulation.write_exact_csv(replication, folder / ("run%04d.csv" % run))
    return simulation.read_exact_csv([folder / ("run%04d.csv" % run) for run in (3, 1, 4, 2)])


def test_read_exact_csv_reads_the_replications_of_r(four_runs, r):
    assert list(four_runs["run"]) == [run for run in (1, 2, 3, 4) for _ in range(8)]
    same(numbers(four_runs), r("simulation_read_exact"))


@pytest.mark.parametrize("kind", ["simulation", "bias"])
def test_summarize_runs(kind, four_runs, r, tmp_path):
    table = simulation.summarize_runs(four_runs, "ZIGP3", simulation.TABLE_COLUMNS[kind])
    assert list(table.columns[:4]) == ["setting", "cutoff", "model", "statistic"]
    assert list(table["statistic"][:2]) == ["mean", "se"] and len(table) == 16
    same(numbers(table), r("simulation_summary_" + kind))
    rbase.write_csv(table, tmp_path / "table.csv")
    assert (tmp_path / "table.csv").read_text().splitlines() \
        == r.lines("simulation_summary_file_" + kind)


def test_summarize_runs_known_values():
    runs = pd.DataFrame({"run": [1.0, 2.0, 3.0], "cutoff": "C1", "model": "P",
                         "x_FDR": [np.nan, np.nan, np.nan], "y": [1.0, 2.0, 6.0]})
    table = simulation.summarize_runs(runs, "S", ["x_FDR", "y"])
    first = table[(table["cutoff"] == "C1") & (table["model"] == "P")]
    assert list(first["y"]) == [3.0, 7.0 ** 0.5]
    assert list(first["x_FDR"]) == [0.0, 0.0]  # an FDR cell without a value is 0
    assert np.isnan(table[table["model"] == "ZIGP"]["y"]).all()  # no replication: NA


def test_setting_tag(r):
    assert simulation.setting_tag("ZIGP1", 5) == "ZIGP1_C5" == r.lines("simulation_setting_tag")[0]


def test_simulation_table_rerun_writes_the_files_of_r(r, tmp_path, capsys):
    """run_replications() and stored_replications() behind a rerun: the
    table and a stored replication, line for line as R writes them."""
    opts = {"rerun": True, "runs": "1:2", "cores": 2}
    paths = {"archive": os.path.join(REPO_ROOT, "archive"), "results": str(tmp_path / "results")}
    table = simulation.simulation_table("Table2", "simulation", "binomial", 0.35, 1, opts, SMALL,
                                        paths)
    same(numbers(table), r("simulation_table"))
    written = tmp_path / "results" / "Table2_simulation_binomial-pi0.35_C1.csv"
    assert written.read_text().splitlines() == r.lines("simulation_table_file")
    run = tmp_path / "results" / "simulation" / "ZIGP4_C1" / "run0002.csv"
    assert run.read_text().splitlines() == r.lines("simulation_run_file")
    assert "6 replications to run" in capsys.readouterr().out

    # finished replications are kept: nothing is left to run the second time
    group = simulation.SIMULATION_SETTINGS.iloc[3:6]

    def run_file(setting, number):
        return str(tmp_path / "results" / "simulation" / simulation.setting_tag(setting, 1)
                   / ("run%04d.csv" % number))

    simulation.run_replications(group, 1, run_file, opts, SMALL)
    assert "0 replications to run" in capsys.readouterr().out
    stored = simulation.stored_replications(group, 1, opts, SMALL, paths)
    assert list(stored) == ["ZIP2", "ZIGP3", "ZIGP4"] and len(stored["ZIGP4"]) == 16


def test_stored_replications_reads_an_archive(tmp_path):
    group = simulation.SIMULATION_SETTINGS.iloc[[1]]
    replication = simulation.simulation_replication(SETTING, 5, 1, SMALL)
    os.makedirs(tmp_path / "archive" / "simulation")
    simulation.write_exact_csv(replication, tmp_path / "run.csv")
    with gzip.open(tmp_path / "archive" / "simulation" / "ZIGP1_C5.csv.gz", "wt") as f:
        f.write((tmp_path / "run.csv").read_text())
    paths = {"archive": str(tmp_path / "archive"), "results": str(tmp_path / "results")}
    stored = simulation.stored_replications(group, 5, {"rerun": False}, SMALL, paths)
    same(numbers(stored["ZIGP1"]), numbers(replication))
    with pytest.raises(SystemExit):  # no archive of this true cutoff
        simulation.stored_replications(group, 1, {"rerun": False}, SMALL, paths)
