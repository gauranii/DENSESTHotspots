# -----------------------------------------------------------------------------
# File:    python/tests/test_setup.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of setup.py: the settings of the paper, the paths, the
#          command line and the messages.
# Usage:   python -m pytest python/tests/test_setup.py
# -----------------------------------------------------------------------------
import os
import re

import pytest

from zinull import setup


def test_settings_of_the_paper():
    assert setup.PAPER == {"seed": 889191, "figure_seed": 889302, "n": 1000, "replications": 1000,
                           "fdr_level": 0.05, "fdr_level_real": 0.10}
    assert os.path.isfile(os.path.join(setup.REPO_ROOT, "flake.nix"))
    assert setup.PATHS["archive"] == os.path.join(setup.REPO_ROOT, "archive")


def test_project_paths(monkeypatch):
    monkeypatch.delenv("ZINULL_OUTPUT", raising=False)
    monkeypatch.delenv("ZINULL_DATA", raising=False)
    paths = setup.project_paths("/repo")
    assert paths == {"data": "/repo/data/realdata.txt", "archive": "/repo/archive",
                     "results": "/repo/results", "figures": "/repo/figures"}
    assert setup.project_paths("/repo", "/out")["figures"] == "/out/figures"
    monkeypatch.setenv("ZINULL_OUTPUT", "/elsewhere")
    monkeypatch.setenv("ZINULL_DATA", "/private/counts.txt")
    paths = setup.project_paths("/repo")
    assert paths["results"] == "/elsewhere/results" and paths["data"] == "/private/counts.txt"
    assert paths["archive"] == "/repo/archive"


def test_require_data(tmp_path, capsys):
    file = tmp_path / "realdata.txt"
    with pytest.raises(SystemExit) as stop:
        setup.require_data(str(file))
    assert stop.value.code != 0 and "ZINULL_DATA" in str(stop.value)
    file.write_text("id\n")
    assert setup.require_data(str(file)) == str(file)


def test_output_file(tmp_path):
    file = setup.output_file(str(tmp_path / "a" / "b"), "t.csv")
    assert file == str(tmp_path / "a" / "b" / "t.csv") and os.path.isdir(tmp_path / "a" / "b")


def test_parse_args():
    defaults = {"runs": "1000", "cores": 4}
    assert setup.parse_args(defaults, []) == {"rerun": False, "runs": "1000", "cores": 4}
    opts = setup.parse_args(defaults, ["--rerun", "--runs=1:2", "--cores=2"])
    assert opts == {"rerun": True, "runs": "1:2", "cores": 2}
    for wrong in (["runs=3"], ["--unknown=1"], ["--runs"], ["--rerun=yes"]):
        with pytest.raises(SystemExit) as stop:
            setup.parse_args(defaults, wrong)
        assert stop.value.code != 0


def test_show_help(tmp_path, capsys):
    script = tmp_path / "script.py"
    script.write_text("# -----\n# File:    script.py\n# Usage:   python script.py\n# -----\nx = 1\n")
    assert setup.show_help([], str(script)) is None
    with pytest.raises(SystemExit) as stop:
        setup.show_help(["--help"], str(script))
    assert stop.value.code == 0
    assert capsys.readouterr().out == "File:    script.py\nUsage:   python script.py\n"


def test_parse_index_set():
    assert setup.parse_index_set("3:7") == [3, 4, 5, 6, 7]
    assert setup.parse_index_set("2,5,9") == [2, 5, 9]
    assert setup.parse_index_set("4") == [4]
    assert setup.parse_index_set("4", count=True) == [1, 2, 3, 4]
    assert setup.parse_index_set("5:3") == [5, 4, 3]


def test_log_message(capsys):
    setup.log_message("%s: %d replications", "Table1", 3)
    assert re.fullmatch(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d - Table1: 3 replications\n",
                        capsys.readouterr().out)
