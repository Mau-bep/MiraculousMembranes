"""Figure style of the manuscript: one module that every plotting script imports.

    import figstyle
    figstyle.use()                                   # once, before the figure is made
    fig, ax = figstyle.figure(cols=1)                # 8.5 cm wide (one column), 3:4 aspect
    ax.plot(x, y)
    ax.set_xlabel(figstyle.label("Energy", "E", "eV"))          # Energy $E$ (eV)
    figstyle.save(fig, "energy.png")                 # checks the rules, writes at 600 dpi

The rules of the journal, and where they are enforced
    width      one column 8.5 cm (COLUMN_CM); cols=1/3, 1/2, 1, 1.5, 2 give that many columns; width_cm=... overrides. The figure is made at its printed size and saved
               WITHOUT bbox_inches="tight" (which would change it), so a font size in points is a size on paper
    text       the smallest capital letter or numeral at least 2 mm tall: Latin Modern digits are 0.644 em tall, so FONT = 9 pt (2.04 mm). Tick labels, axis labels, legend,
               annotations: all the same size
    lines      linewidth at least 0.5 pt (LW_MIN); the default is 1 pt
    markers    diameter at least 1 mm (2.84 pt, MARKER_MIN_PT); the default is 4 pt; filled, no small open symbols; no hatching, no shading
    check(fig) walks a figure and reports whatever breaks these rules (save() calls it)

Text is LaTeX (Latin Modern, the set-up that was in Plots.py) when latex, dvipng and ghostscript are installed, mathtext with Computer Modern otherwise; use(usetex=False) forces the second.
"""
import glob
import math
import os
import re
import shutil
import warnings

import matplotlib as mpl

# ---- the rules ---------------------------------------------------------------------------------------------------------------------------------------------
COLUMN_CM = 8.5                       # width of one manuscript column (3 3/8 in)
MM_PT = 72.0 / 25.4                   # points per mm
MIN_CAP_MM = 2.0                      # smallest capital letter / numeral
DIGIT_EM = 0.644                      # height of a Latin Modern (Computer Modern) digit in em; capitals are 0.683
FONT = 9.0                            # pt: 9 * 0.644 = 5.80 pt = 2.04 mm
MIN_FONT = MIN_CAP_MM * MM_PT / DIGIT_EM          # 8.80 pt
LW = 1.0                              # default curve width, pt
LW_MIN = 0.5                          # pt (0.18 mm)
AXES_LW = 0.6                         # axes, ticks, legend frame
MARKER = 4.0                          # default marker diameter, pt (1.41 mm)
MARKER_MIN_PT = 1.0 * MM_PT           # 2.84 pt
SCATTER_S = MARKER ** 2               # scatter's s is the area in pt^2 (a square of side = diameter)
DPI = 600
FIG_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Results", "Paper_figures")     # where the scripts dump the figure iterations

_PREAMBLE = R"""
\usepackage{lmodern}
\usepackage{amsmath}
\usepackage{graphicx}
\DeclareMathOperator{\arcsinh}{arcsinh}
\DeclareMathOperator{\km}{k_\mathrm{m}}
\DeclareMathOperator{\fbi}{f_\mathrm{bi}}
\DeclareMathOperator{\eps}{\epsilon_{\mathrm{mc}}}
\DeclareMathOperator{\epscrit}{\epsilon_{\mathrm{mc}}^*}
\DeclareMathOperator{\uf}{u_{\mathrm{f}}}
\DeclareMathOperator{\kBT}{k_\mathrm{B}T}
"""[1:]


def width_cm(cols=1, width=None):
    """Width in cm of a figure `cols` columns wide (1/3, 1/2, 1, 1.5, 2, any number); width (cm) overrides."""
    return float(width) if width is not None else float(cols) * COLUMN_CM


def latex_available():
    return all(shutil.which(p) for p in ("latex", "dvipng", "gs"))


def use(usetex=None, font=FONT):
    """Load the style into matplotlib's rcParams (idempotent). usetex None: LaTeX text when it is installed. font in pt, not below MIN_FONT."""
    if font < MIN_FONT - 1e-9:
        warnings.warn("figstyle: font %.2f pt is under the %.2f pt that gives a %g mm numeral" % (font, MIN_FONT, MIN_CAP_MM))
    if usetex is None:
        usetex = latex_available()
    rc = {
        "text.usetex": bool(usetex),
        "font.family": "serif",
        "font.size": font,
        "axes.labelsize": font, "axes.titlesize": font, "figure.titlesize": font, "figure.labelsize": font,
        "xtick.labelsize": font, "ytick.labelsize": font, "legend.fontsize": font, "legend.title_fontsize": font,
        "axes.linewidth": AXES_LW,
        "xtick.major.width": AXES_LW, "ytick.major.width": AXES_LW, "xtick.minor.width": AXES_LW * 0.8, "ytick.minor.width": AXES_LW * 0.8,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.minor.size": 1.75, "ytick.minor.size": 1.75,
        "xtick.direction": "in", "ytick.direction": "in",
        "lines.linewidth": LW, "lines.markersize": MARKER, "lines.markeredgewidth": 0.0,
        "patch.linewidth": AXES_LW, "grid.linewidth": 0.4,
        "scatter.edgecolors": "none",
        "legend.frameon": False, "legend.borderaxespad": 0.4, "legend.handlelength": 1.6, "legend.labelspacing": 0.3,
        "axes.grid": False,
        "figure.constrained_layout.use": True, "figure.constrained_layout.w_pad": 2.0 / 72, "figure.constrained_layout.h_pad": 2.0 / 72,
        "savefig.dpi": DPI, "savefig.bbox": "standard", "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if usetex:
        rc["text.latex.preamble"] = _PREAMBLE
    else:
        rc.update({"font.serif": ["cmr10"], "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True, "axes.unicode_minus": False})
    mpl.rcParams.update(rc)
    return bool(usetex)


def figure(cols=1, aspect=0.75, height_cm=None, width=None, nrows=1, ncols=1, **kw):
    """A figure of the printed size: cols columns wide (or width cm) and aspect = height / width (or height_cm). Returns fig, axes as plt.subplots does.
    The extra keywords go to plt.subplots (sharex, gridspec_kw, subplot_kw, ...)."""
    import matplotlib.pyplot as plt
    w = width_cm(cols, width) / 2.54
    h = (height_cm / 2.54) if height_cm is not None else w * aspect
    return plt.subplots(nrows, ncols, figsize=(w, h), **kw)


def label(name=None, symbol=None, unit=None, exp=None):
    """An axis label in the journal style: the quantity in words, its symbol in italics, the unit in parentheses.
        label("Energy", "E", "eV")                -> Energy $E$ (eV)
        label("Energy", "E", "eV", exp=7)         -> Energy $E$ ($10^{7}$ eV)           the factor in the label, so the ticks are 0.2, 0.4, ...
        label(symbol=r"\\Delta R/R")              -> $\\Delta R/R$                      a dimensionless quantity: no unit, no parentheses
        label(symbol=r"\\Delta R/R", exp=4)       -> $10^{4}\\,\\Delta R/R$              the factor in front of the symbol
        label("Number of particles")              -> Number of particles                a text label (dimensionless): roman, no unit
    A unit with maths in it is written with $...$ by the caller, e.g. unit=r"$k_\\mathrm{B}T$"."""
    if symbol is None:
        return name or ""
    if unit is None:
        sym = symbol if exp is None else r"10^{%d}\,%s" % (exp, symbol)
        return "%s $%s$" % (name, sym) if name else "$%s$" % sym
    head = "%s $%s$" % (name, symbol) if name else "$%s$" % symbol
    return head + " (%s)" % (unit if exp is None else "$10^{%d}$ %s" % (exp, unit))


def log_label(symbol, unit=None, base=10):
    """The label of a logarithmic axis, the base always given (never a bare log): log_label("E", "MeV") -> $\\log_{10}|E\\,(\\mathrm{MeV})|$"""
    inner = symbol if unit is None else r"%s\,(\mathrm{%s})" % (symbol, unit)
    return r"$\log_{%s}|%s|$" % ("e" if base == "e" else base, inner) if base != "e" else r"$\ln|%s|$" % inner


def panel_label(ax, text, x=0.0, y=1.0, dx=2.0, dy=-2.0):
    """A panel letter ((a), (b), ...) at the top left corner of the axes, offset (dx, dy) pt, in the same size as the rest of the text."""
    return ax.annotate(text, xy=(x, y), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points", ha="left", va="top")


def _is_open(marker, face, edge):
    return marker not in (None, "", "None", " ") and str(face).lower() == "none"


def check(fig, tol=1e-6):
    """The rules of the journal applied to a drawn figure: a list of strings, one per violation (empty: all fine).
    Checked: the size of every text, the width of every line, the diameter of every marker and every scatter point, open markers, hatching.
    Not checked: the size of sub/superscripts inside a LaTeX string (TeX sets them at 70 %), nor an image."""
    from matplotlib.collections import PathCollection
    from matplotlib.lines import Line2D
    from matplotlib.text import Text
    out = []
    fig.canvas.draw()
    for t in fig.findobj(Text):
        s = t.get_text()
        if s.strip() and t.get_visible() and t.get_fontsize() < MIN_FONT - tol:
            out.append("text %r is %.2f pt (minimum %.2f)" % (s[:30], t.get_fontsize(), MIN_FONT))
    furniture = set()                 # tick marks and grid lines belong to the axes, not to the data
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            furniture.update(id(l) for l in axis.get_ticklines(minor=False) + axis.get_ticklines(minor=True) + axis.get_gridlines())
    for ln in fig.findobj(Line2D):
        if not ln.get_visible() or id(ln) in furniture:
            continue
        if ln.get_linestyle() not in ("None", "", " ") and ln.get_linewidth() < LW_MIN - tol:
            out.append("line %r is %.2f pt wide (minimum %.2f)" % (ln.get_label()[:30], ln.get_linewidth(), LW_MIN))
        if ln.get_marker() not in (None, "", "None", " ", ",") and len(ln.get_xdata()) > 0:
            if ln.get_markersize() < MARKER_MIN_PT - tol:
                out.append("marker of %r is %.2f pt across (minimum %.2f)" % (ln.get_label()[:30], ln.get_markersize(), MARKER_MIN_PT))
            if _is_open(ln.get_marker(), ln.get_markerfacecolor(), ln.get_markeredgecolor()):
                out.append("marker of %r is open (use filled symbols)" % ln.get_label()[:30])
    for c in fig.findobj(PathCollection):
        if not c.get_visible():
            continue
        sizes = c.get_sizes()
        if len(sizes) and min(sizes) < MARKER_MIN_PT ** 2 - tol:
            out.append("scatter %r: s = %.2f pt^2, a diameter of %.2f pt (minimum %.2f)" % (c.get_label()[:30], min(sizes), math.sqrt(min(sizes)), MARKER_MIN_PT))
    for p in fig.findobj(mpl.patches.Patch):
        if p.get_hatch():
            out.append("hatching on %s (avoid hatching and shading)" % type(p).__name__)
    return out


def save(fig, path, dpi=None, check_rules=True, **kw):
    """Write the figure at its printed size (no bbox_inches="tight") and report violations of the rules. PNG at 600 dpi by default; a path ending .pdf is vector."""
    if kw.pop("bbox_inches", None) not in (None, "standard"):
        warnings.warn("figstyle.save: bbox_inches is ignored, it would change the printed size")
    if check_rules:
        for msg in check(fig):
            warnings.warn("figstyle: " + msg)
    fig.savefig(path, dpi=dpi or DPI, **kw)
    return path


def save_figure(fig, name, folder=None, formats=("png", "eps"), version=True, **kw):
    """Write the figure to the figure folder (FIG_DIR = Results/Paper_figures) in every format of `formats` (png and eps by default), at its printed size, and report violations of the rules.
    version=True numbers the iterations: name_v01.png/.eps, name_v02.png/.eps, ... (the next number after the highest in the folder), so an earlier iteration is never overwritten;
    version=False overwrites name.png/.eps. Returns the list of the files written."""
    folder = folder or FIG_DIR
    os.makedirs(folder, exist_ok=True)
    stem = name
    if version:
        used = [int(m.group(1)) for f in glob.glob(os.path.join(folder, name + "_v*.*")) for m in [re.search(r"_v(\d+)\.[a-z]+$", f)] if m]
        stem = "%s_v%02d" % (name, max(used, default=0) + 1)
    if kw.pop("check_rules", True):
        for msg in check(fig):
            warnings.warn("figstyle: " + msg)
    out = []
    for fmt in formats:
        path = os.path.join(folder, "%s.%s" % (stem, fmt))
        fig.savefig(path, dpi=DPI, **kw)
        out.append(path)
    print("wrote " + ", ".join(out))
    return out
