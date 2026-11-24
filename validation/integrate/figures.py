"""Figures of docs/9_validation_time_integration.md (matplotlib, static PNG in docs/images/).

Made by report.py from the stored results (and, for the curves and the example, from tac_benchmarks).
Reference palette of the dataviz method (validated: categorical slots 1-4 on the light surface, with
direct labels because slots 3-4 are below 3:1 contrast)."""

from pathlib import Path

import numpy as np

import tac_benchmarks as B

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # categorical slots 1-4 (fixed order)

IMG_DIR = Path(__file__).resolve().parents[2] / "docs" / "images" / "validation_integrate"
FIGURES = {
    "curves": "curves.png",
    "bias": "bias.png",
    "coverage": "coverage.png",
    "example": "example.png",
}


def _style():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "axes.titlesize": 10, "axes.labelsize": 9, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.6, "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
        "ytick.labelcolor": INK2, "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "legend.frameon": False, "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.5, "svg.hashsalt": "radot", "path.simplify": True,
    })
    return plt


def _save(fig, key):
    IMG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(IMG_DIR / FIGURES[key], dpi=150, bbox_inches="tight", metadata={"Software": None})


def fig_curves(plt):
    """The three exact curves, with the times of each schedule."""
    t = np.linspace(0, 260, 1000)
    fig, axes = plt.subplots(1, 3, figsize=(9, 2.6), sharex=True)
    for ax, (name, s) in zip(axes, B.SHAPES.items()):
        ax.plot(t, s.curve(t), color=SERIES[0])
        ax.set_title(f"{name}: {s.description}", fontsize=8.5, loc="left")
        ax.set_xlabel("time after injection (h)")
        top = float(np.max(s.curve(t)))
        for i, (sc, times) in enumerate(B.SCHEDULES.items()):
            y = -0.08 * top * (i + 1)
            ax.plot(times, [y] * len(times), "|", color=INK2, markersize=7, mew=1.2)
            ax.text(262, y, sc, va="center", fontsize=7, color=INK2)
        ax.set_ylim(-0.32 * top, 1.08 * top)
        ax.set_xlim(0, 260)
        ax.set_yticks([v for v in ax.get_yticks() if 0 <= v <= 1.08 * top])
        ax.spines["left"].set_bounds(0, 1.08 * top)
    axes[0].set_ylabel("activity (MBq)")
    fig.tight_layout()
    _save(fig, "curves")
    plt.close(fig)


def fig_bias(plt, res, noise="0.1"):
    """Mean relative error ± SD per method: small multiples shape (rows) x schedule (columns)."""
    methods = list(B.METHODS)
    shapes, schedules = list(res["results"]), list(res["settings"]["schedules"])
    fig, axes = plt.subplots(len(shapes), len(schedules), figsize=(9, 7.5), sharey=True)
    y = np.arange(len(methods))[::-1]
    for i, sh in enumerate(shapes):
        for j, sc in enumerate(schedules):
            ax = axes[i, j]
            ax.axvline(0, color=INK2, linewidth=0.9)
            for k, m in enumerate(methods):
                r = res["results"][sh][sc][noise][m]
                if not r["applicable"] or r["bias_pct"] is None:
                    ax.text(-57, y[k], "n/a", fontsize=7, color=MUTED, va="center", ha="left")
                    continue
                ax.errorbar(r["bias_pct"], y[k], xerr=r["sd_pct"], fmt="o", color=SERIES[0], markersize=4,
                            elinewidth=1.2, capsize=0)
            ax.set_xlim(-60, 80)
            ax.grid(axis="y", visible=False)
            if i == 0:
                ax.set_title(f"{sc} ({len(res['settings']['schedules'][sc])} points)")
            if j == 0:
                ax.set_yticks(y, methods)
                ax.text(-0.95, 0.5, sh, transform=ax.transAxes, rotation=90, va="center", fontsize=10,
                        color=INK, fontweight="bold")
            if i == len(shapes) - 1:
                ax.set_xlabel("TIA error (%)")
    fig.suptitle(f"Mean error ± SD of the TIA, {100 * float(noise):g} % noise", x=0.02, ha="left", fontsize=10,
                 color=INK)
    fig.tight_layout()
    _save(fig, "bias")
    plt.close(fig)


def fig_coverage(plt, res):
    """Coverage of ±1 SE for every configuration with a standard error, by method; 68 % expected."""
    methods = [m for m in B.METHODS if any(
        "coverage_68" in res["results"][sh][sc][no][m]
        for sh in res["results"] for sc in res["results"][sh] for no in res["results"][sh][sc])]
    fig, ax = plt.subplots(figsize=(7, 3.0))
    y = np.arange(len(methods))[::-1]
    lo, hi = B.COVERAGE_68_RANGE
    ax.axvspan(lo, hi, color=GRID, zorder=0)
    ax.axvline(0.68, color=INK2, linewidth=0.9, zorder=1)
    ax.text(0.68, len(methods) - 0.35, "68 % expected (band: CI acceptance)", fontsize=7.5, color=INK2,
            ha="center")
    rng = np.random.default_rng(0)
    shapes = list(res["results"])
    for c, sh in enumerate(shapes):  # 3 curves: categorical slots 1-3 (valid for all pairs)
        for k, m in enumerate(methods):
            vals = [res["results"][sh][sc][no][m]["coverage_68"] for sc in res["results"][sh]
                    for no in res["results"][sh][sc] if "coverage_68" in res["results"][sh][sc][no][m]]
            ax.plot(vals, y[k] + rng.uniform(-0.2, 0.2, len(vals)), "o", color=SERIES[c], markersize=5,
                    markeredgecolor=SURFACE, markeredgewidth=0.8, label=sh if k == 0 else None)
    ax.legend(loc="lower left", title="curve", title_fontsize=8)
    ax.set_yticks(y, methods)
    ax.set_xlim(0, 1)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{100 * v:.0f} %"))
    ax.set_xlabel("fraction of realisations with |TIA − exact| ≤ SE (one dot per curve × schedule × noise)")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    _save(fig, "coverage")
    plt.close(fig)


def fig_example(plt, shape="lesion", schedule="4tp", noise=0.10, index=0):
    """One noisy realisation and the curves of three methods, against the exact curve."""
    from radot.integrate import TAC_MODELS, evaluate

    s = B.SHAPES[shape]
    times = np.asarray(B.SCHEDULES[schedule])
    values = B.realisations(s, schedule, noise, index + 1)[index]
    t = np.linspace(0, 300, 1200)
    fig, ax = plt.subplots(figsize=(7, 3.4))
    exact = s.tia_mbq_h()
    ax.plot(t, s.curve(t), color=INK, linewidth=1.5, linestyle=(0, (4, 3)))
    series = []
    label_at = {"f2 fsd": (75, -0.22), "f4a fsd": (55, 0.18), "trapezoid ramp/phys": (215, 0.14)}
    for i, label in enumerate(("f2 fsd", "f4a fsd")):
        r = B.run_one(B.METHODS[label], times, values, noise)
        fit = r["fits"][r["model"]]
        model = TAC_MODELS[r["model"]]
        p = np.array([fit["params"][n] for n in model.param_names])
        series.append((label, evaluate(model, p, t, B.LAMBDA_P), r["tia_mbq_h"], SERIES[i]))
    tr = B.run_one(B.METHODS["trapezoid ramp/phys"], times, values, noise)
    a_n = values[-1]
    trap = np.interp(t, np.r_[0.0, times], np.r_[0.0, values])
    trap[t > times[-1]] = a_n * np.exp(-B.LAMBDA_P * (t[t > times[-1]] - times[-1]))
    series.append(("trapezoid ramp/phys", trap, tr["tia_mbq_h"], SERIES[2]))
    for label, curve, tia, color in series:
        ax.plot(t, curve, color=color, label=f"{label}: {100 * (tia / exact - 1):+.0f} %")
        x, dy = label_at[label]
        ax.text(x, float(np.interp(x, t, curve)) + dy, label, color=INK2, fontsize=7.5, ha="center")
    ax.plot(times, values, "o", color=INK, markersize=6, markeredgecolor=SURFACE, markeredgewidth=1.0,
            label=f"measured ({100 * noise:g} % noise)", zorder=5)
    ax.plot([], [], color=INK, linestyle=(0, (4, 3)), label="exact curve")
    ax.set_xlim(0, 300)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("time after injection (h)")
    ax.set_ylabel("activity (MBq)")
    ax.set_title(f"{shape}, {schedule}: one noisy realisation (legend: TIA error)", loc="left")
    ax.legend(loc="upper right")
    fig.tight_layout()
    _save(fig, "example")
    plt.close(fig)


def make_all(res):
    plt = _style()
    fig_curves(plt)
    fig_bias(plt, res)
    fig_coverage(plt, res)
    fig_example(plt)
