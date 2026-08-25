"""
Figure 2 — BIRDMAn use cases, rendered in full.

    a  Vaginal microbiome (Ravel et al. 2011): BIRDMAn posterior differentials
       for the N most BV-associated and N most health-associated features.
       Bar colour ramps with |differential| within each direction. Only features
       whose 95% HDI excludes zero are eligible; ranking on the posterior mean
       alone promotes features whose means are large only because they are
       barely estimated.

    b  Dual-course ciprofloxacin (Dethlefsen & Relman 2011): per-subject sample
       log-ratio trajectories for the FirstCp and SecondCp contrasts (top row),
       and the corresponding posterior derivative with a 5-95% ribbon (bottom).

Both panels run from data committed to this repository. Panel b is drawn by
scripts/relman_abx/3.01-plot_figure2b.py so it keeps the styling of the original
published figure; if its inputs under results/relman_abx/ are ever absent, a
labelled placeholder is drawn so the composite still renders.

HDI unpacking follows https://github.com/ahdilmore/MARS_Birdman
(birdman/zebra_birdman_analysis.ipynb): split on the comma, strip the
parentheses, flag credible when the interval excludes zero, then convert the
bounds to offsets from the mean for errorbar().

Usage:
    python scripts/figure2.py [-n 10] [--include-non-credible] [--no-ci] [--seed 0]
"""

import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Patch

REPO = Path(__file__).resolve().parents[1]

# panel b is drawn by the relman_abx script so the composite keeps the styling
# of the original published figure rather than restyling it here
_FIG2B = SourceFileLoader(
    "fig2b", str(REPO / "scripts" / "relman_abx" / "3.01-plot_figure2b.py")
).load_module()

# authored at final print width so production does not scale the type down
FIG_W_MM, FIG_H_MM = 180.0, 114.0          # OUP double column
MM = 1 / 25.4
BASE_PT = 6.5                              # "large" -> 7.8, "x-large" -> 9.4

INK, MUTED, RULE, GRID = "#222222", "#666666", "#BBBBBB", "#DDDDDD"

# one set of chrome weights for both panels, so a and b read as one figure.
# OUP's floor is 0.25 pt, so nothing here may go below it.
SPINE_LW = 0.7
TICK_LW = 0.6
TICK_LEN = 2.5
GRID_LW = 0.4
ERR_LW = 0.7

# notebooks/paper.mplstyle, minus the savefig/dpi keys -- those are global and
# would override this figure's explicit GridSpec margins. grid.linewidth is
# raised from grahman's 0.1, which is below the printable minimum. he sized
# panel b with relative keywords, so font.size rescales it as a whole.
PANEL_B_STYLE = {
    "figure.facecolor": "white",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "font.size": BASE_PT,
    "axes.labelsize": "large",
    "axes.titlesize": "x-large",
    "axes.titlelocation": "left",
    "axes.axisbelow": True,
    "grid.linewidth": GRID_LW,
    "grid.color": GRID,
    # absolute, not grahman's "x-small": relative to a 6.5pt base that would
    # render the subject legend at 4.5pt
    "legend.fontsize": BASE_PT - 0.5,
    "legend.title_fontsize": BASE_PT,
    "lines.linewidth": 1.0,
    "xtick.labelsize": BASE_PT,
    "ytick.labelsize": BASE_PT,
    # shared chrome -- must match panel_a below
    "axes.edgecolor": RULE,
    "axes.linewidth": SPINE_LW,
    "xtick.color": RULE,
    "ytick.color": RULE,
    "xtick.labelcolor": INK,
    "ytick.labelcolor": INK,
    "xtick.major.width": TICK_LW,
    "ytick.major.width": TICK_LW,
    "xtick.major.size": TICK_LEN,
    "ytick.major.size": TICK_LEN,
}
DIFFERENTIALS = REPO / "data" / "qadabra" / "outputs" / "birdman-differentials.tsv"
RELMAN = REPO / "results" / "relman_abx"
RELMAN_TBL = REPO / "data" / "relman_abx" / "processed" / "processed_tbl.biom"
RELMAN_MD = REPO / "data" / "relman_abx" / "processed" / "processed_md.tsv"
OUTDIR = REPO / "figures"

COVARIATE = "C(study_condition, Treatment('healthy'))[T.bacterial_vaginosis]"
LEVELS = ["preCp", "FirstCp", "FirstWPC", "Interim", "SecondCp", "SecondWPC", "PostCp"]
LEVEL_DIFFS = [f"{LEVELS[i]}_vs_{LEVELS[i - 1]}" for i in range(1, len(LEVELS))]
CONTRASTS = ["FirstCp_vs_preCp", "SecondCp_vs_Interim"]
N_LOGRATIO = 40

# Panel a: two sequential ramps about a pale midpoint (a diverging scheme keyed
# to magnitude). Panel b subjects: categorical, validated colourblind-safe --
# chroma >= 0.1, adjacent CVD dE 10.5 (protan), normal-vision dE 26.9.
BV_CMAP, HEALTH_CMAP = plt.get_cmap("Reds"), plt.get_cmap("Blues")
RAMP_LO, RAMP_HI = 0.22, 0.92
BV_KEY, HEALTH_KEY = BV_CMAP(0.75), HEALTH_CMAP(0.55)


# ---------------------------------------------------------------------------
# Panel a — vaginal differentials
# ---------------------------------------------------------------------------
def unpack_hdi(df, covariate):
    """'(lo, hi)' -> errorbar offsets from the mean, plus a credible flag."""
    hdi_col, mean_col = f"{covariate}_hdi", f"{covariate}_mean"
    bounds = df[hdi_col].astype(str).str.split(",", expand=True)
    df["lower"] = bounds[0].str.strip().str[1:].astype(float)
    df["upper"] = bounds[1].str.strip().str[:-1].astype(float)
    df["credible"] = np.where((df["lower"] > 0) | (df["upper"] < 0), "yes", "no")
    df["upper"] = df["upper"] - df[mean_col]
    df["lower"] = df[mean_col] - df["lower"]
    return df


def ramp_colors(means):
    """Shade each bar by |differential|, scaled within its own direction."""
    out = []
    pos, neg = means[means > 0], means[means < 0]
    pmax = pos.max() if len(pos) else 1.0
    nmax = abs(neg.min()) if len(neg) else 1.0
    for m in means:
        frac = abs(m) / (pmax if m > 0 else nmax)
        shade = RAMP_LO + (RAMP_HI - RAMP_LO) * frac
        out.append(BV_CMAP(shade) if m > 0 else HEALTH_CMAP(shade))
    return out


def panel_a(ax, n, credible_only, show_ci):
    df = (pd.read_csv(DIFFERENTIALS, sep="\t")
            .rename(columns={"feature id": "Feature"})
            .set_index("Feature"))
    df = unpack_hdi(df, COVARIATE)
    mean_col = f"{COVARIATE}_mean"
    if credible_only:
        df = df[df["credible"] == "yes"]

    df = df.sort_values(mean_col)
    sel = df if len(df) <= 2 * n else pd.concat([df.head(n), df.tail(n)])

    ys = np.arange(len(sel))          # ascending value, so most positive on top
    means = sel[mean_col].to_numpy()

    ax.barh(ys, means, color=ramp_colors(means), height=0.74, zorder=2)
    if show_ci:
        ax.errorbar(means, ys, xerr=sel[["lower", "upper"]].to_numpy().T,
                    fmt="none", ecolor=MUTED, elinewidth=ERR_LW, capsize=1.5,
                    capthick=ERR_LW, alpha=0.9, zorder=3)
    ax.axvline(0, color=RULE, linewidth=SPINE_LW, zorder=4)

    ax.set_yticks(ys)
    ax.set_yticklabels([f.replace("_", " ") for f in sel.index], fontsize=BASE_PT)
    ax.set_ylim(-0.7, len(sel) - 0.3)
    # the differential is the posterior coefficient; matches panel b's beta bar
    ax.set_xlabel(r"BIRDMAn $\beta$", fontsize=BASE_PT + 1.5, color=INK,
                  labelpad=4)
    ax.tick_params(axis="both", labelsize=BASE_PT, length=TICK_LEN,
                   width=TICK_LW, color=RULE, labelcolor=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(RULE)
        ax.spines[side].set_linewidth(SPINE_LW)

    # Legend below the panel, horizontal — as in the reference layout.
    ax.legend(handles=[Patch(facecolor=BV_KEY, label="BV-associated"),
                       Patch(facecolor=HEALTH_KEY, label="Healthy-associated")],
              loc="upper center", bbox_to_anchor=(0.5, -0.125), ncol=2,
              fontsize=BASE_PT, frameon=False, handlelength=1.2,
              handleheight=0.9, columnspacing=1.4, borderpad=0.0)


# ---------------------------------------------------------------------------
# Panel b — antibiotic log-ratio dynamics
# ---------------------------------------------------------------------------
def relman_available():
    """panel b needs the committed model outputs; report whether they are here."""
    beta = RELMAN / "birdman_results.beta_var.tsv"
    quants = [RELMAN / f"{c}_derivative_quantiles.tsv" for c in CONTRASTS]
    return (beta.exists() and RELMAN_TBL.exists() and RELMAN_MD.exists()
            and all(q.exists() for q in quants))


def panel_b_placeholder(axes):
    for ax in axes:
        ax.set_axis_off()
    ax = axes[0]
    ax.add_patch(FancyBboxPatch(
        (0.02, -1.18), 2.14, 2.16, boxstyle="round,pad=0,rounding_size=0.02",
        transform=ax.transAxes, facecolor="#F5F6F8", edgecolor=RULE,
        linewidth=0.7, linestyle=(0, (3, 2)), clip_on=False, zorder=1))
    ax.text(1.09, 0.20, "Antibiotic log-ratio dynamics", transform=ax.transAxes,
            ha="center", va="center", fontsize=10, fontweight="bold", color=INK)
    ax.text(1.09, -0.08,
            "awaiting BIRDMAn model outputs\n\n"
            "results/relman_abx/birdman_results.beta_var.tsv\n"
            "results/relman_abx/beta_var.nc  (or *_derivative_quantiles.tsv)",
            transform=ax.transAxes, ha="center", va="center", fontsize=7.5,
            color=MUTED, linespacing=1.7)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=10, help="features per side in panel a")
    ap.add_argument("--include-non-credible", action="store_true",
                    help="do not require the HDI to exclude zero in panel a")
    ap.add_argument("--no-ci", action="store_true",
                    help="omit HDI error bars from panel a (reference styling)")
    ap.add_argument("--seed", type=int, default=None,
                    help="seed panel b's bootstrap CI bands for reproducibility")
    args = ap.parse_args()

    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
    # wspace here is the gap between panel a and panel b; panel b's own internal
    # spacing is set by 3.01. it needs to clear panel b's math ylabel.
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.22, 1.0, 1.0],
                  hspace=0.14, wspace=0.62, left=0.170, right=0.988,
                  top=0.90, bottom=0.17)
    ax_a = fig.add_subplot(gs[:, 0])
    panel_a(ax_a, args.n, credible_only=not args.include_non_credible,
            show_ci=not args.no_ci)

    if relman_available():
        with plt.rc_context(PANEL_B_STYLE):
            _FIG2B.plot_figure_2b(fig=fig, subplot_spec=gs[:, 1:],
                                  seed=args.seed)
    else:
        print(f"[panel b] missing model outputs under {RELMAN.relative_to(REPO)}"
              " - drawing placeholder")
        panel_b_placeholder([fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]),
                             fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])])

    # Panel labels in figure coordinates so long tick labels cannot push them
    # off-canvas.
    for x, label in ((0.010, "a"), (0.455, "b")):
        fig.text(x, 0.935, label, fontsize=BASE_PT + 4, fontweight="bold",
                 va="bottom", ha="left", color=INK)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        out = OUTDIR / f"figure2.{ext}"
        fig.savefig(out, dpi=450 if ext == "png" else None, facecolor="white")
        print("wrote", out.relative_to(REPO))
    plt.close(fig)


if __name__ == "__main__":
    main()
