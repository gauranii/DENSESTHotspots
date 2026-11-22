# -----------------------------------------------------------------------------
# File:    python/tests/test_real_data.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of real_data.py on the counts of two simulated "domains"
#          (python/tests/data/small_domains.txt, in the layout of
#          data/realdata.txt): reading, null fits and the two tables, against
#          R, double for double.
# Usage:   python -m pytest python/tests/test_real_data.py
# -----------------------------------------------------------------------------
import pytest

from conftest import numbers, same
from zinull import real_data


def test_read_domains(r, data_file, tmp_path):
    domains = real_data.read_domains(data_file)
    assert list(domains) == ["domA", "domB"]
    sizes = [len(domains["domA"]), len(domains["domB"]), domains["domA"].sum(),
             domains["domB"].sum()]
    same(sizes, r("real_domains"))  # the NA padding of domB is left out
    assert domains["domA"].dtype.kind == "i"
    with pytest.raises(SystemExit):
        real_data.read_domains(tmp_path / "missing.txt")


def test_fit_domain(r, data_file):
    x = real_data.read_domains(data_file)["domA"]
    fits = real_data.fit_domain(x)
    assert (fits["D"] <= x.max()).all()  # the bound D is at most the largest count
    same(numbers(fits), r("real_fit_domain"))


def test_protein_rejections(r, data_file):
    table = real_data.protein_rejections(data_file, 0.10)
    assert list(table.columns) == ["domain", "cutoff", "procedure", "ZIGP", "ZIP", "GP", "P"]
    assert list(table["procedure"][:3]) == ["one_stage", "two_stage", "storey"]
    assert len(table) == 12
    same(numbers(table), r("real_rejections"))


def test_protein_estimates(r, data_file):
    table = real_data.protein_estimates(data_file)
    assert list(table.columns) == ["domain", "cutoff", "model", "eta", "lambda", "theta", "pi0",
                                   "C", "D"]
    assert len(table) == 16 and list(table["domain"].unique()) == ["domA", "domB"]
    same(numbers(table), r("real_estimates"))
