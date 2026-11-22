# -----------------------------------------------------------------------------
# File:    python/tests/test_plots.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Biometrics 74:458-471.
# Purpose: Tests of plots.py: the bars, points and lines of each figure are
#          the plotted values, in R's order and colours, and the files are
#          written under R's names.
# Usage:   python -m pytest python/tests/test_plots.py
# -----------------------------------------------------------------------------
import os

import numpy as np
import pandas as pd
from matplotlib.colors import to_hex

from zinull import plots
from zinull.figure_data import bar_heights, null_fit_series
from zinull.real_data import fit_domain, read_domains

NULL = np.array([30.0, 12, 5, 1])
NON_NULL = np.array([0.0, 0, 2, 4, 3, 1])
SETTINGS = {"n": 60, "figure_seed": 889302}


def bars(figure):
    """(left, width, height, fill, border) of every bar, in drawing order."""
    return [(p.get_x(), p.get_width(), p.get_height(), to_hex(p.get_facecolor()),
             to_hex(p.get_edgecolor())) for p in figure.axes[0].patches]


def test_draw_non_null_first():
    drawn = bars(plots.draw_non_null_first(NULL, NON_NULL))
    # the non-null counts in black first, then the null counts in gray with gray borders
    assert [bar[2] for bar in drawn] == list(NON_NULL) + list(NULL)
    assert [bar[0] for bar in drawn] == [0, 1, 2, 3, 4, 5, 0, 1, 2, 3]
    assert {bar[3:] for bar in drawn[:6]} == {("#000000", "#000000")}
    assert {bar[3:] for bar in drawn[6:]} == {("#bebebe", "#bebebe")}
    assert all(bar[1] == 1 for bar in drawn)
    assert plots.draw_non_null_first(NULL, NON_NULL).axes[0].get_ylim() == (0, 50)


def test_draw_null_first():
    figure = plots.draw_null_first(NULL, NON_NULL)
    drawn = bars(figure)
    assert [bar[2] for bar in drawn] == list(NULL) + list(NON_NULL)
    assert {bar[3:] for bar in drawn[:4]} == {("#bebebe", "#000000")}
    assert {bar[3:] for bar in drawn[4:]} == {("#000000", "#000000")}
    assert figure.axes[0].get_ylim() == (0, 50)
    assert figure.get_size_inches().tolist() == [7, 7]


def test_draw_protein_histogram():
    drawn = bars(plots.draw_protein_histogram(bar_heights([0, 0, 1, 3])))
    assert [bar[2] for bar in drawn] == [2, 1, 0, 1]
    # barplot() leaves a space of 0.2 before every bar of width 1
    np.testing.assert_allclose([bar[0] for bar in drawn], [0.2, 1.4, 2.6, 3.8])
    assert {bar[3:] for bar in drawn} == {("#bebebe", "#000000")}


def test_draw_null_fit():
    series = pd.DataFrame({"fEstimation": ["zip", "zip", "fHat", "fHat", "fHat", "zigp", "zigp"],
                           "C": [0, 1, 0, 1, 2, 0, 1],
                           "values": [0.5, 0.2, 0.45, 0.25, 0.1, 0.4, 0.3]})
    lines = plots.draw_null_fit(series).axes[0].lines
    assert [line.get_label() for line in lines] == ["fHat", "zigp", "zip"]
    assert [line.get_marker() for line in lines] == ["o", "^", "s"]
    assert lines[0].get_xdata().tolist() == [0, 1, 2]
    assert lines[0].get_ydata().tolist() == [0.45, 0.25, 0.1]
    assert lines[2].get_ydata().tolist() == [0.5, 0.2]
    assert all(line.is_dashed() for line in lines)


def test_save_figure(tmp_path, capsys):
    plots.save_figure("bars.pdf", str(tmp_path / "figures"),
                      lambda: plots.draw_null_first(NULL, NON_NULL))
    content = (tmp_path / "figures" / "bars.pdf").read_bytes()
    assert content.startswith(b"%PDF") and b"/MediaBox [ 0 0 504 504 ]" in content  # 7 x 7 in
    assert "wrote figures/bars.pdf" in capsys.readouterr().out


def test_save_simulation_histograms(tmp_path):
    paths = {"figures": str(tmp_path)}
    plots.save_simulation_histograms("Fig1", SETTINGS, paths)
    plots.save_simulation_histograms("FigS3", SETTINGS, paths)
    assert sorted(os.listdir(tmp_path)) == [
        "Fig1A_simulation_histogram_ZIP1-geometric.pdf",
        "Fig1B_simulation_histogram_ZIGP1-geometric.pdf",
        "Fig1C_simulation_histogram_ZIGP2-geometric.pdf",
        "Fig1D_simulation_histogram_ZIP2-binomial.pdf",
        "Fig1E_simulation_histogram_ZIGP3-binomial.pdf",
        "Fig1F_simulation_histogram_ZIGP4-binomial.pdf",
        "FigS3A_simulation_histogram_MZIGP1-geometric.pdf",
        "FigS3B_simulation_histogram_MZIGP2-geometric.pdf",
        "FigS3C_simulation_histogram_MZIP1-geometric.pdf"]


def test_save_protein_histograms(tmp_path, data_file):
    paths = {"figures": str(tmp_path), "data": data_file}
    plots.save_protein_histograms("Fig2", ["domB", "domA"], paths)
    assert sorted(os.listdir(tmp_path)) == ["Fig2A_protein_histogram_domB.pdf",
                                            "Fig2B_protein_histogram_domA.pdf"]


def test_save_null_fit_plots(tmp_path, data_file):
    paths = {"figures": str(tmp_path), "data": data_file}
    plots.save_null_fit_plots("FigS7", paths, domains={"FigS7": "domA"})
    assert sorted(os.listdir(tmp_path)) == ["FigS7A_protein_null-fit_domA_C1.pdf",
                                            "FigS7B_protein_null-fit_domA_C2.pdf"]
    x = read_domains(data_file)["domA"]
    series = null_fit_series(x, fit_domain(x), "C2")
    lines = plots.draw_null_fit(series).axes[0].lines
    assert sum(len(line.get_xdata()) for line in lines) == len(series)
