"""
Figure 2 — BIRDMAn use cases, rendered in full.

    a  Vaginal microbiome (Ravel et al. 2011): posterior differentials with 95%
       HDI error bars, for the 10 most BV-associated and 10 most
       health-associated features. Only features whose HDI excludes zero are
       eligible — ranking on the posterior mean alone promotes features whose
       means are large only because they are barely estimated.

    b  Dual-course ciprofloxacin (Dethlefsen & Relman 2011): per-subject sample
       log-ratio trajectories with 95% bootstrap CI bands (top row), and the
       posterior derivative with a 5-95% ribbon (bottom row), for the FirstCp
       and SecondCp contrasts.

Both panels run from committed artifacts, so this works from a fresh clone with
no cluster and no model refit:
    data/qadabra/outputs/birdman-differentials.tsv
    results/relman_abx/birdman_results.beta_var.tsv
    results/relman_abx/{contrast}_derivative_quantiles.tsv
    data/relman_abx/processed/{processed_tbl.biom,processed_md.tsv}

The figure is authored at final print width (180 mm), so the sizes in the file
are the sizes on the page.

HDI unpacking follows https://github.com/ahdilmore/MARS_Birdman: split on the
comma, strip the parentheses, flag credible when the interval excludes zero,
then convert the bounds to offsets from the mean for errorbar().

Usage:
    python scripts/figure2.py
"""

from pathlib import Path

import biom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

REPO = Path(__file__).resolve().parents[1]
DIFFERENTIALS = REPO / "data" / "qadabra" / "outputs" / "birdman-differentials.tsv"
RELMAN = REPO / "results" / "relman_abx"
RELMAN_TBL = REPO / "data" / "relman_abx" / "processed" / "processed_tbl.biom"
RELMAN_MD = REPO / "data" / "relman_abx" / "processed" / "processed_md.tsv"
OUTDIR = REPO / "figures"

COVARIATE = "C(study_condition, Treatment('healthy'))[T.bacterial_vaginosis]"
LEVELS = ["preCp", "FirstCp", "FirstWPC", "Interim", "SecondCp", "SecondWPC", "PostCp"]
LEVEL_DIFFS = [f"{LEVELS[i]}_vs_{LEVELS[i - 1]}" for i in range(1, len(LEVELS))]
CONTRASTS = ["FirstCp_vs_preCp", "SecondCp_vs_Interim"]
SUBJECTS = ["494.D", "494.E", "494.F"]
N_PER_SIDE = 10      # features shown per direction in panel a
N_LOGRATIO = 40      # features in each half of the panel b log-ratio
N_BOOT, SEED = 1000, 0    # bootstrap CI for the panel b bands; seeded, so the
                          # bands are identical run to run

FIG_W_MM, FIG_H_MM = 180.0, 114.0          # OUP double column
MM = 1 / 25.4
BASE_PT = 6.5

# Panel a: two sequential ramps about a pale midpoint. Panel b subjects:
# categorical, validated colourblind-safe (protan dE 9.2, normal dE 22.9).
BV_CMAP, HEALTH_CMAP = plt.get_cmap("Reds"), plt.get_cmap("Blues")
RAMP = (0.22, 0.92)
SUBJ_COLORS = ["#0173b2", "#de8f05", "#029e73"]
INK, MUTED, RULE, GRID = "#222222", "#666666", "#BBBBBB", "#DDDDDD"

# OUP's line-weight floor is 0.25 pt.
SPINE_LW, TICK_LW, TICK_LEN, GRID_LW, ERR_LW = 0.7, 0.6, 2.5, 0.4, 0.7

# beta_bar_bottom < 0 for every contrast, so the betas cannot go inside the log.
_RATIO = (r"\ln\left(\frac{\sum\ \mathrm{Top\ 40\ OTUs}}"
          r"{\sum \mathrm{Bottom\ 40\ OTUs}}\right)")
TOP_LABEL = f"${_RATIO}$"
BOT_LABEL = rf"$\Delta\ {_RATIO}$"


def style_axis(ax, grid=False):
    """Shared chrome, so both panels read as one figure."""
    if grid:
        ax.grid(True, color=GRID, linewidth=GRID_LW)
        ax.set_axisbelow(True)
    ax.tick_params(labelsize=BASE_PT, length=TICK_LEN, width=TICK_LW,
                   color=RULE, labelcolor=INK)
    for spine in ax.spines.values():
        spine.set_color(RULE)
        spine.set_linewidth(SPINE_LW)


# ---------------------------------------------------------------------------
# Panel a — vaginal differentials
# ---------------------------------------------------------------------------
def panel_a(ax):
    mean_col, hdi_col = f"{COVARIATE}_mean", f"{COVARIATE}_hdi"
    df = (pd.read_csv(DIFFERENTIALS, sep="\t")
            .rename(columns={"feature id": "Feature"})
            .set_index("Feature"))

    lo, hi = (df[hdi_col].astype(str).str.strip("()")
              .str.split(",", expand=True).astype(float).to_numpy().T)
    df["lower"] = df[mean_col] - lo          # errorbar wants offsets, not bounds
    df["upper"] = hi - df[mean_col]
    df = df[(lo > 0) | (hi < 0)]             # credible: HDI excludes zero

    df = df.sort_values(mean_col)
    sel = pd.concat([df.head(N_PER_SIDE), df.tail(N_PER_SIDE)])
    means = sel[mean_col].to_numpy()
    ys = np.arange(len(sel))                 # ascending, so most positive on top

    scale = {1: means.max(), -1: abs(means.min())}
    colors = [(BV_CMAP if m > 0 else HEALTH_CMAP)(
                  RAMP[0] + (RAMP[1] - RAMP[0]) * abs(m) / scale[np.sign(m)])
              for m in means]

    ax.barh(ys, means, color=colors, height=0.74, zorder=2)
    ax.errorbar(means, ys, xerr=sel[["lower", "upper"]].to_numpy().T, fmt="none",
                ecolor=MUTED, elinewidth=ERR_LW, capsize=1.5, capthick=ERR_LW,
                alpha=0.9, zorder=3)
    ax.axvline(0, color=RULE, linewidth=SPINE_LW, zorder=4)

    ax.set_yticks(ys)
    ax.set_yticklabels([f.replace("_", " ") for f in sel.index], fontsize=BASE_PT)
    ax.set_ylim(-0.7, len(sel) - 0.3)
    ax.set_xlabel(r"BIRDMAn $\beta$", fontsize=BASE_PT + 1.5, color=INK, labelpad=4)
    style_axis(ax)
    for side, spine in ax.spines.items():
        spine.set_visible(side in ("left", "bottom"))

    ax.legend(handles=[Patch(facecolor=BV_CMAP(0.75), label="BV-associated"),
                       Patch(facecolor=HEALTH_CMAP(0.55),
                             label="Healthy-associated")],
              loc="upper center", bbox_to_anchor=(0.5, -0.125), ncol=2,
              fontsize=BASE_PT, frameon=False, handlelength=1.2,
              handleheight=0.9, columnspacing=1.4, borderpad=0.0)


# ---------------------------------------------------------------------------
# Panel b — antibiotic log-ratio dynamics
# ---------------------------------------------------------------------------
def load_relman():
    summ = pd.read_table(RELMAN / "birdman_results.beta_var.tsv", index_col=0)
    summ.index = summ.index.astype(str)
    # patsy names backward-difference columns after the first level of each
    # pair, so the rename is positional -- column order is load-bearing
    summ.columns = (["Intercept_mean"] + [f"{d}_mean" for d in LEVEL_DIFFS]
                    + ["Intercept_std"] + [f"{d}_std" for d in LEVEL_DIFFS])
    cent = (summ[[f"{d}_mean" for d in LEVEL_DIFFS]]
            .apply(lambda x: x - x.mean(), axis=0))

    table = biom.load_table(str(RELMAN_TBL)).to_dataframe(dense=True).T
    md = pd.read_table(RELMAN_MD, index_col=0)
    md["antibiotic"] = pd.Categorical(md["antibiotic"], categories=LEVELS,
                                      ordered=True)

    log_ratios = {}
    for contrast in CONTRASTS:
        ranked = cent[f"{contrast}_mean"].sort_values(ascending=False)
        num = table[ranked.head(N_LOGRATIO).index].sum(axis=1) + 1
        den = table[ranked.tail(N_LOGRATIO).index].sum(axis=1) + 1
        log_ratios[contrast] = np.log(num / den)

    derivs = {}
    for contrast in CONTRASTS:
        d = pd.read_table(RELMAN / f"{contrast}_derivative_quantiles.tsv")
        d["covariate"] = pd.Categorical(d["covariate"], categories=LEVEL_DIFFS,
                                        ordered=True)
        derivs[contrast] = d.sort_values("covariate")

    return pd.DataFrame(log_ratios).join(md), derivs


def boot_ci(values, rng):
    """95% bootstrap CI of the mean — the band seaborn draws, but seeded."""
    if len(values) < 2:
        return np.nan, np.nan
    draws = rng.choice(values, size=(N_BOOT, len(values)), replace=True)
    return np.percentile(draws.mean(axis=1), [2.5, 97.5])


def panel_b(fig, gs):
    lr_df, derivs = load_relman()
    rng = np.random.default_rng(SEED)
    colors = dict(zip(SUBJECTS, SUBJ_COLORS))
    x_of = {lvl: i for i, lvl in enumerate(LEVELS)}

    lr1 = fig.add_subplot(gs[0, 0])
    lr2 = fig.add_subplot(gs[0, 1], sharey=lr1)
    for ax, contrast in zip((lr1, lr2), CONTRASTS):
        for subj, g in lr_df.groupby("host_subject_id", observed=True):
            pts = g.groupby("antibiotic", observed=True)[contrast]
            xs = [x_of[lvl] for lvl, _ in pts]
            mean = [v.mean() for _, v in pts]
            ci = np.array([boot_ci(v.to_numpy(), rng) for _, v in pts])
            ax.fill_between(xs, ci[:, 0], ci[:, 1], color=colors[subj],
                            alpha=0.20, linewidth=0, zorder=2)
            ax.plot(xs, mean, marker="o", ms=3, lw=1.0, color=colors[subj],
                    markeredgecolor="white", markeredgewidth=0.4, zorder=3)
        ax.set_title(contrast.replace("_vs_", " vs. "), fontsize=BASE_PT + 1.5)
        ax.set_xticks(range(len(LEVELS)))
        ax.set_xticklabels([])
        ax.set_xlim(-0.25, 6.25)
        style_axis(ax, grid=True)
    lr1.set_ylabel(TOP_LABEL, fontsize=BASE_PT + 2.5)
    lr2.tick_params("y", width=0, labelleft=False)
    lr1.legend(handles=[Patch(facecolor=colors[s], label=s, linewidth=0.4,
                              edgecolor=INK) for s in SUBJECTS],
               title="Subject", title_fontsize=BASE_PT, fontsize=BASE_PT - 0.5,
               loc="upper right", frameon=True, edgecolor=RULE,
               labelspacing=0.3, handlelength=1.2)

    d1 = fig.add_subplot(gs[1, 0])
    d2 = fig.add_subplot(gs[1, 1], sharey=d1)
    for ax, contrast in zip((d1, d2), CONTRASTS):
        d = derivs[contrast]
        lo, mid, hi = (d[d["quantile"] == q]["beta_var"].to_numpy()
                       for q in (0.05, 0.5, 0.95))
        xs = np.arange(len(mid))
        ax.fill_between(xs, lo, hi, color="gray", alpha=0.25, linewidth=0, zorder=2)
        for bound in (lo, hi):
            ax.plot(xs, bound, color="gray", ls="--", lw=0.8, zorder=3)
        ax.plot(xs, mid, marker="o", ms=3, lw=1.0, color="gray",
                markeredgecolor=INK, markeredgewidth=0.4, zorder=4)
        ax.set_xticks(np.arange(-1, 6))
        ax.set_xticklabels(LEVELS, rotation=45, ha="right", fontsize=BASE_PT)
        ax.set_xlim(-1.25, 5.25)
        style_axis(ax, grid=True)
    d1.set_ylabel(BOT_LABEL, fontsize=BASE_PT + 2.5)
    d2.tick_params("y", width=0, labelleft=False)


# ---------------------------------------------------------------------------
def main():
    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
    # wspace has to clear panel b's math ylabel
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.22, 1.0, 1.0],
                  hspace=0.18, wspace=0.62, left=0.170, right=0.988,
                  top=0.90, bottom=0.17)
    panel_a(fig.add_subplot(gs[:, 0]))
    panel_b(fig, gs[:, 1:].subgridspec(2, 2, wspace=0.08, hspace=0.18))

    for x, label in ((0.010, "a"), (0.455, "b")):
        fig.text(x, 0.935, label, fontsize=BASE_PT + 4, fontweight="bold",
                 va="bottom", ha="left", color=INK)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext, dpi in (("pdf", None), ("png", 450)):
        out = OUTDIR / f"figure2.{ext}"
        fig.savefig(out, dpi=dpi, facecolor="white")
        print("wrote", out.relative_to(REPO))
    plt.close(fig)


if __name__ == "__main__":
    main()
