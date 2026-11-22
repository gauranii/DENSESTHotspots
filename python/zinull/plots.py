# -----------------------------------------------------------------------------
# File:    python/zinull/plots.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. Figures 1 and 2, Web Figures S1 to S12.
# Purpose: Histograms of counts: null and non-null counts of a simulated data
#          set drawn over each other, and the counts of a protein domain;
#          observed against fitted null probabilities of a protein domain.
#          The Python counterpart of R/plots.R (matplotlib): the same bars,
#          points and lines, on pages laid out as R's barplot() and ggplot2
#          lay them out.
# Inputs:  None.
# Outputs: figures/<name>.pdf (written by save_figure()).
# Usage:   from zinull.setup import PAPER, PATHS
#          from zinull.plots import save_simulation_histograms
#          save_simulation_histograms("Fig1", PAPER, PATHS)
# Source:  R/plots.R.
# -----------------------------------------------------------------------------
import logging
import math

import matplotlib
import numpy as np
from matplotlib.backends.backend_pdf import FigureCanvasPdf
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

from .cutoff import CUTOFFS
from .figure_data import (HISTOGRAM_SETTINGS, NULL_FIT_DOMAINS, bar_heights, null_fit_series,
                          protein_counts, simulate_histogram_data)
from .real_data import fit_domain, read_domains
from .setup import log_message, output_file

PAGE = 7.0     # width and height of a figure, inches (the default of R's pdf())
GRAY = "#BEBEBE"  # R's "gray"
LINE = 0.75    # line width of R's base graphics, points
# Plot region of R's base graphics on the page (margins of 5.1, 4.1, 4.1 and
# 2.1 lines of 0.2 inches), as fractions of the page: left, bottom, width, height.
BASE_REGION = (0.82 / PAGE, 1.02 / PAGE, (PAGE - 0.82 - 0.42) / PAGE, (PAGE - 1.02 - 0.82) / PAGE)
# Panel of the ggplot2 figures on the page.
GGPLOT_PANEL = (40.98 / 504, 31.91 / 504, 381.21 / 504, 466.61 / 504)
# The page content is written as text (not compressed), so that the marks of
# a figure can be read back; the fonts are the PDF core fonts, as in R.
STYLE = {"pdf.compression": 0, "pdf.use14corefonts": True, "font.family": "Helvetica"}
LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


# ----- Histograms (R's barplot()) -----

def _base_axes(xlim, ylim, cex):
    """A page with the plot region of R's base graphics: the horizontal
    range widened by 4 percent at each end, the vertical one as given, and
    an axis on the left only."""
    figure = Figure(figsize=(PAGE, PAGE))
    axes = figure.add_axes(BASE_REGION)
    pad = 0.04 * (xlim[1] - xlim[0])
    axes.set_xlim(xlim[0] - pad, xlim[1] + pad)
    axes.set_ylim(*ylim)
    for side in ("top", "right", "bottom"):
        axes.spines[side].set_visible(False)
    axes.spines["left"].set_linewidth(LINE)
    axes.yaxis.set_major_locator(MaxNLocator(nbins=7, steps=[1, 2, 5, 10]))
    ticks = [tick for tick in axes.get_yticks() if 0 <= tick <= ylim[1]]
    axes.set_yticks(ticks)
    axes.spines["left"].set_bounds(ticks[0], ticks[-1])
    axes.tick_params(axis="y", width=LINE, length=7.2, labelsize=12 * cex, labelrotation=90)
    for label in axes.get_yticklabels():
        label.set_verticalalignment("center")
    axes.tick_params(axis="x", length=0, labelsize=12 * cex, pad=14)
    axes.set_xticks([])
    return figure, axes


def _bars(axes, heights, space, fill, border, cex, names=True):
    """Bars as `barplot(heights, space = space, col = fill, border = border)`
    places them (width 1), with the counts under the bars that have room."""
    heights = np.asarray(heights, dtype=np.float64)
    left = space + (1 + space) * np.arange(len(heights))
    axes.bar(left, heights, 1.0, align="edge", color=fill, edgecolor=border, linewidth=LINE,
             clip_on=False)
    if names:
        per_bar = (1 + space) * axes.get_window_extent().width / np.diff(axes.get_xlim())[0]
        widest = 0.556 * 12 * cex * len(str(len(heights) - 1)) * 100 / 72  # pixels at 100 dpi
        step = max(1, math.ceil(1.6 * widest / per_bar))
        shown = np.arange(0, len(heights), step)
        axes.set_xticks(left[shown] + 0.5, [str(i) for i in shown])
    return left


def draw_non_null_first(null, non_null):
    """Histogram with the non-null counts (black) under the null counts.

    Parameters: null, non_null (bar heights of the null and the non-null
    counts).
    Returns: the figure (matplotlib Figure).
    """
    null, non_null = np.asarray(null, dtype=np.float64), np.asarray(non_null, dtype=np.float64)
    figure, axes = _base_axes((0, len(non_null)), (0, np.max(null) + 20), cex=0.8)
    _bars(axes, non_null, 0, "black", "black", cex=0.8)
    _bars(axes, null, 0, GRAY, GRAY, cex=0.8, names=False)
    return figure


def draw_null_first(null, non_null):
    """Histogram with the null counts (gray) under the non-null counts (black).

    Parameters: null, non_null (bar heights of the null and the non-null
    counts).
    Returns: the figure (matplotlib Figure).
    """
    null, non_null = np.asarray(null, dtype=np.float64), np.asarray(non_null, dtype=np.float64)
    largest_count = max(len(null), len(non_null)) - 1
    figure, axes = _base_axes((0, largest_count), (0, np.max(null) + 20), cex=0.7)
    _bars(axes, null, 0, GRAY, "black", cex=0.7, names=len(null) >= len(non_null))
    _bars(axes, non_null, 0, "black", "black", cex=0.7, names=len(non_null) > len(null))
    return figure


def draw_protein_histogram(heights):
    """Histogram of the counts of a protein domain.

    Parameters: heights (bar heights, in the order of the counts 0, 1, ...).
    Returns: the figure (matplotlib Figure).
    """
    heights = np.asarray(heights, dtype=np.float64)
    top = np.max(heights)
    figure, axes = _base_axes((0.2, 1.2 * len(heights)), (-0.01 * top, top), cex=0.7)
    _bars(axes, heights, 0.2, GRAY, "black", cex=0.7)
    return figure


# ----- Null fits (ggplot2) -----

def draw_null_fit(series):
    """Observed and fitted null probabilities against the count.

    Parameters: series (output of `null_fit_series()`).
    Returns: the figure (matplotlib Figure): one dashed line with points
    per series (circles, triangles, squares, as ggplot2 gives them).
    """
    figure = Figure(figsize=(PAGE, PAGE))
    axes = figure.add_axes(GGPLOT_PANEL)
    names = sorted(series["fEstimation"].unique())
    sizes = {"o": 5.0, "^": 6.7, "s": 5.0}  # as ggplot2 draws its shapes at size 2
    for name, marker in zip(names, ("o", "^", "s")):
        part = series[series["fEstimation"] == name].sort_values("C", kind="stable")
        axes.plot(part["C"], part["values"], marker=marker, markersize=sizes[marker], color="black",
                  markeredgewidth=0, linewidth=0.32, dashes=(4, 4), label=name, clip_on=False)
    axes.margins(0.05)
    axes.grid(True, color="#EBEBEB", linewidth=1.07)
    axes.set_axisbelow(True)
    axes.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    axes.yaxis.set_major_locator(MaxNLocator(nbins=5))
    for spine in axes.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(1.07)
    axes.tick_params(colors="#4D4D4D", labelsize=9, length=2.75, width=1.07)
    axes.set_xlabel("C", fontsize=11)
    axes.set_ylabel("estimated probabilities", fontsize=11)
    axes.yaxis.set_label_coords(-0.068, 0.5)
    axes.legend(title="fEstimation", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False,
                fontsize=9, title_fontsize=11, numpoints=1, alignment="left")
    return figure


# ----- Files -----

def save_figure(name, figures_dir, draw):
    """Save a figure as pdf (7 x 7 in, the default of R's `pdf()`).

    Parameters: name (file name), figures_dir (folder of the figures), draw
    (function without arguments that returns the figure).
    Returns: None.
    """
    fonts = logging.getLogger("matplotlib.font_manager")
    level = fonts.level
    fonts.setLevel(logging.ERROR)  # no notes about the weights of the core fonts
    try:
        with matplotlib.rc_context(STYLE):
            FigureCanvasPdf(draw()).print_pdf(output_file(figures_dir, name))
    finally:
        fonts.setLevel(level)
    log_message("wrote figures/%s", name)


def save_simulation_histograms(figure, settings, paths, panels=HISTOGRAM_SETTINGS):
    """Simulate and save the histograms of one figure.

    The files are named <figure><panel>_simulation_histogram_<model>-<signal>.pdf.

    Parameters: figure (figure of the paper, for example "FigS3"), settings
    (settings of the paper: `n`, `figure_seed`), paths (project paths:
    `figures`), panels (table of the panels, as `HISTOGRAM_SETTINGS`).
    Returns: None.
    """
    panels = panels[panels["figure"] == figure]
    for _, panel in panels.iterrows():
        data = simulate_histogram_data(panel, settings["n"], settings["figure_seed"])
        null = bar_heights(data["count"][data["non_null"] == 0])
        non_null = bar_heights(data["count"][data["non_null"] == 1])
        draw = {"non_null": draw_non_null_first, "null": draw_null_first}[panel["order"]]
        name = "%s%s_simulation_histogram_%s-%s.pdf" % (figure, panel["panel"], panel["model"],
                                                        panel["signal"])
        save_figure(name, paths["figures"], lambda: draw(null, non_null))


def save_protein_histograms(figure, domains, paths):
    """Save the histograms of the protein domains.

    The files are named <figure><panel>_protein_histogram_<domain>.pdf.

    Parameters: figure (figure of the paper, "Fig2"), domains (names of the
    domains, in the order of the panels), paths (project paths: `data`,
    `figures`).
    Returns: None.
    """
    for i, domain in enumerate(domains):
        heights = bar_heights(protein_counts(domain, paths["data"]))
        name = "%s%s_protein_histogram_%s.pdf" % (figure, LETTERS[i], domain)
        save_figure(name, paths["figures"], lambda: draw_protein_histogram(heights))


def save_null_fit_plots(figure, paths, domains=NULL_FIT_DOMAINS):
    """Save the null-fit plots of one protein domain, one per cutoff criterion.

    The files are named <figure><panel>_protein_null-fit_<domain>_<criterion>.pdf
    (panel A: C1, panel B: C2).

    Parameters: figure (figure of the paper, for example "FigS7"), paths
    (project paths: `data`, `figures`), domains (domain of each figure, as
    `NULL_FIT_DOMAINS`).
    Returns: None.
    """
    domain = domains[figure]
    x = read_domains(paths["data"])[domain]
    fits = fit_domain(x)
    for i, criterion in enumerate(CUTOFFS):
        name = "%s%s_protein_null-fit_%s_%s.pdf" % (figure, LETTERS[i], domain, criterion)
        series = null_fit_series(x, fits, criterion)
        save_figure(name, paths["figures"], lambda: draw_null_fit(series))
