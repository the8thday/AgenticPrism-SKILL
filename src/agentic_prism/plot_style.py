"""Shared figure styling. Appearance only: no function here changes data, limits or fitted values."""
import logging
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cycler import cycler

# CJK fallback fonts often lack a bold face; matplotlib then substitutes weight 600 and logs each lookup.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# Okabe-Ito order with black first: colour-blind safe and legible in greyscale print.
# Yellow (#F0E442) is left out because it has too little contrast on white.
OKABE_ITO = ["#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F"]
STANDARD = ["#355db2", "#b5651d", "#6a4c93", "#2a9d8f", "#c1121f", "#606c38", "#3a86ff", "#8d99ae"]

THEMES = {
    "prism_like": {"color": "#000000", "marker": "o", "grid": False, "axes_width": 1.2, "line_width": 1.6,
                   "palette": OKABE_ITO, "markers": ["o", "s", "^", "D", "v", "p", "h", "<"],
                   "spine_offset_pt": 5, "bold_labels": True, "residual_marker": "open"},
    "standard": {"color": "#355db2", "marker": "o", "grid": True, "axes_width": .8, "line_width": 1.5,
                 "palette": STANDARD, "markers": ["o"] * 8,
                 "spine_offset_pt": 0, "bold_labels": False, "residual_marker": "filled"},
}

MUTED = "#8a8f95"      # reference lines, excluded points, extrapolations
EXCLUDED = "#b74d46"   # explicitly excluded observations


def rc(token, font, cjk=None):
    """rcParams for one theme; use as ``with plt.rc_context(rc(token, font, cjk))``."""
    weight = "bold" if token["bold_labels"] else "normal"
    width = token["axes_width"]
    return {
        "font.family": [font] + ([cjk] if cjk else []), "font.size": 9,
        "axes.labelsize": 10, "axes.labelweight": weight, "axes.titlesize": 9.5, "axes.titleweight": weight,
        "axes.titlelocation": "left", "axes.linewidth": width, "axes.edgecolor": "black",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": token["grid"], "grid.alpha": .18, "grid.linewidth": .6, "axes.axisbelow": True,
        "axes.prop_cycle": cycler(color=token["palette"]),
        "xtick.direction": "out", "ytick.direction": "out", "xtick.color": "black", "ytick.color": "black",
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "xtick.major.width": width, "ytick.major.width": width, "xtick.minor.width": width * .75, "ytick.minor.width": width * .75,
        "xtick.major.size": 4.5, "ytick.major.size": 4.5, "xtick.minor.size": 2.5, "ytick.minor.size": 2.5,
        "lines.linewidth": token["line_width"], "lines.markersize": 5, "errorbar.capsize": 3,
        "legend.frameon": False, "legend.fontsize": 8, "legend.handlelength": 1.6, "legend.borderaxespad": .3,
        "pdf.fonttype": 42, "svg.fonttype": "path", "axes.unicode_minus": False,
        "savefig.facecolor": "white", "figure.facecolor": "white",
    }


def color(token, index):
    palette = token["palette"]
    return palette[index % len(palette)]


def marker(token, index):
    markers = token["markers"]
    return markers[index % len(markers)]


def point_style(token, index=0, size=5):
    """Keyword arguments for ``ax.plot(..., linestyle='none', **point_style(...))`` data markers."""
    c = color(token, index)
    # Prism draws filled symbols with a thin black outline so overlapping replicates stay distinguishable.
    edge = "black" if token["bold_labels"] and c != "#000000" else c
    return {"marker": marker(token, index), "ms": size, "mfc": c, "mec": edge, "mew": .6, "linestyle": "none"}


def residual_style(token, index=0):
    c = color(token, index)
    if token["residual_marker"] == "open":
        return {"marker": "o", "ms": 3.6, "mfc": "white", "mec": c, "mew": .9, "linestyle": "none"}
    return {"marker": "o", "ms": 3.6, "mfc": c, "mec": c, "alpha": .75, "linestyle": "none"}


def display_label(name):
    """Readable axis text for declared column names; identifiers with spaces are kept verbatim."""
    name = str(name)
    return name.replace("_", " ") if " " not in name else name


def _is_plot_axes(ax):
    # Heatmaps, colorbars and hidden axes keep their default frame.
    return ax.axison and not ax.images and not hasattr(ax, "_colorbar") and ax.get_label() != "<colorbar>"


def finish(fig, token):
    """Detach visible left/bottom spines outward as Prism does, so axes do not meet at the origin."""
    offset = token["spine_offset_pt"]
    for ax in fig.axes:
        if not _is_plot_axes(ax):
            continue
        for side in ("left", "bottom"):
            if offset and ax.spines[side].get_visible():
                ax.spines[side].set_position(("outward", offset))
    return fig


def save(fig, stem, token, dpi, formats=("svg", "pdf", "png")):
    """Finish and write ``stem.<ext>`` for each format, then close the figure."""
    finish(fig, token)
    for ext in formats:
        fig.savefig(f"{stem}.{ext}", dpi=dpi, facecolor="white")
    plt.close(fig)


def significance_label(p):
    """GraphPad-style summary of an already saved, already adjusted p value."""
    if p is None:
        return None
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < .05 else "ns"


def significance_brackets(ax, positions, pairs, top, span):
    """Draw brackets for (i, j, p) tuples above ``top``; returns the new upper y limit.

    Positions are x coordinates of groups; ``span`` is the data range used to size gaps.
    Only saved adjusted p values are drawn; no test is performed here.
    """
    step, tick = span * .09, span * .02
    y = top + step * .6
    for i, j, p in sorted(pairs, key=lambda t: (abs(t[1] - t[0]), min(t[0], t[1]))):
        text = significance_label(p)
        if text is None:
            continue
        x1, x2 = positions[i], positions[j]
        ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick], color="black", lw=.9, clip_on=False)
        ax.text((x1 + x2) / 2, y + (0 if text != "ns" else tick * .4), text, ha="center", va="bottom",
                fontsize=8.5 if text != "ns" else 7.5)
        y += step
    return y + step * .3


def swarm_offsets(values, y_span, width=.32):
    """Deterministic beeswarm: each point takes the nearest free slot left or right of the centre."""
    import numpy as np
    values = np.asarray(values, float)
    offsets = np.zeros(len(values))
    if len(values) < 2 or not np.isfinite(y_span) or y_span <= 0:
        return offsets
    near, step = y_span * .035, width / 6
    placed = []
    for i in np.argsort(values, kind="stable"):
        for k in range(0, 13):
            candidate = (k + 1) // 2 * step * (1 if k % 2 else -1)
            if abs(candidate) > width / 2:
                candidate = ((i % 5) - 2) * step / 2
                break
            if all(abs(values[i] - v) >= near or abs(candidate - o) >= step * .95 for v, o in placed):
                break
        offsets[i] = candidate
        placed.append((values[i], candidate))
    return offsets


def dot_plot(ax, token, groups, center="mean", spread="sd"):
    """Prism-style scatter dot plot. ``groups`` is a list of (label, values); returns x positions.

    Every independent value is drawn. ``center`` is "mean" or "median"; ``spread`` is "sd" or None.
    """
    import numpy as np
    finite = np.concatenate([np.asarray(v, float) for _, v in groups]) if groups else np.array([])
    finite = finite[np.isfinite(finite)]
    span = float(finite.max() - finite.min()) if len(finite) else 1.
    positions = np.arange(len(groups), dtype=float)
    prism = token["bold_labels"]
    for i, (label, values) in enumerate(groups):
        v = np.asarray(values, float)
        x = positions[i] + swarm_offsets(v, span)
        if prism:
            ax.plot(x, v, zorder=3, **point_style(token, i, 5.5))
        else:
            ax.plot(x, v, marker="o", ms=5, linestyle="none", color=token["color"], alpha=.8, zorder=3)
        if not len(v):
            continue
        middle = float(np.mean(v)) if center == "mean" else float(np.median(v))
        ax.hlines(middle, positions[i] - .24, positions[i] + .24, color="black", lw=1.6 if prism else 1.4, zorder=4)
        if spread == "sd" and len(v) > 1:
            ax.errorbar(positions[i], middle, yerr=float(np.std(v, ddof=1)), fmt="none", ecolor="black",
                        elinewidth=1.1, capsize=6, capthick=1.1, zorder=4)
    ax.set_xticks(positions, [str(label) for label, _ in groups])
    if len(groups) > 4 or max((len(str(label)) for label, _ in groups), default=0) > 10:
        for t in ax.get_xticklabels():
            t.set_rotation(35)
            t.set_ha("right")
            t.set_rotation_mode("anchor")
    ax.set_xlim(-.6, len(groups) - .4)
    ax.tick_params(axis="x", length=0 if not prism else 4.5)
    return positions
