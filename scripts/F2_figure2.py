"""
Figure 2 — BIRDMAn use cases.

    a  Vaginal microbiome (Ravel et al. 2011): posterior differentials with 95%
       HDI error bars, for the 10 most BV-associated and 10 most
       health-associated features. Only features whose HDI excludes zero are
       eligible — ranking on the posterior mean alone promotes features whose
       means are large only because they are barely estimated.

    b  Dual-course ciprofloxacin (Dethlefsen & Relman 2011): per-subject sample
       log-ratio trajectories (top row) and the posterior derivative with a
       5-95% ribbon (bottom row), for the FirstCp and SecondCp contrasts.

Panel a runs from committed data. Panel b also needs the BIRDMAn model outputs
under results/relman_abx/, which are not committed; a labelled placeholder is
drawn when they are absent.

HDI unpacking follows https://github.com/ahdilmore/MARS_Birdman: split on the
comma, strip the parentheses, flag credible when the interval excludes zero,
then convert the bounds to offsets from the mean for errorbar().

Usage:
    python scripts/F2_figure2.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

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
N_PER_SIDE = 10      # features shown per direction in panel a
N_LOGRATIO = 40      # features in each half of the panel b log-ratio

# Panel a: two sequential ramps about a pale midpoint. Panel b subjects:
# categorical, validated colourblind-safe (protan dE 10.5, normal dE 26.9).
BV_CMAP, HEALTH_CMAP = plt.get_cmap("Reds"), plt.get_cmap("Blues")
RAMP = (0.22, 0.92)
SUBJ_COLORS = ["#0077BB", "#EE7733", "#009988"]
INK, MUTED, RULE, GRID = "#222222", "#666666", "#BBBBBB", "#DDDDDD"


# ---------------------------------------------------------------------------
# Panel a
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

    # Shade each bar by |differential|, scaled within its own direction.
    scale = {1: means.max(), -1: abs(means.min())}
    colors = [(BV_CMAP if m > 0 else HEALTH_CMAP)(
                  RAMP[0] + (RAMP[1] - RAMP[0]) * abs(m) / scale[np.sign(m)])
              for m in means]

    ax.barh(ys, means, color=colors, height=0.74, zorder=2)
    ax.errorbar(means, ys, xerr=sel[["lower", "upper"]].to_numpy().T, fmt="none",
                ecolor=INK, elinewidth=0.6, capsize=1.5, capthick=0.6,
                alpha=0.8, zorder=3)
    ax.axvline(0, color=INK, linewidth=0.8, zorder=4)

    ax.set_yticks(ys)
    ax.set_yticklabels([f.replace("_", " ") for f in sel.index], fontsize=8)
    ax.set_ylim(-0.7, len(sel) - 0.3)
    ax.set_xlabel("BIRDMAn Differential", fontsize=9.5, color=INK, labelpad=6)
    ax.tick_params(labelsize=8, length=3, color=RULE)
    for side, spine in ax.spines.items():
        spine.set_visible(side in ("left", "bottom"))
        spine.set_color(RULE)

    ax.legend(handles=[Patch(facecolor=BV_CMAP(0.75), label="BV-associated"),
                       Patch(facecolor=HEALTH_CMAP(0.55), label="Healthy-associated")],
              loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2,
              fontsize=9, frameon=False, handlelength=1.4, columnspacing=2.0)


# ---------------------------------------------------------------------------
# Panel b
# ---------------------------------------------------------------------------
def load_relman():
    """Pieces panel b needs, or None if the model outputs are absent."""
    beta_tsv = RELMAN / "birdman_results.beta_var.tsv"
    if not (beta_tsv.exists() and RELMAN_TBL.exists() and RELMAN_MD.exists()):
        return None

    import biom

    summ = pd.read_table(beta_tsv, sep="\t", index_col=0)
    summ.index = summ.index.astype(str)
    # Columns are Intercept, 6 means, Intercept_std, 6 stds — positional, per
    # scripts/relman_abx/2.01-summarize_birdman.py.
    summ.columns = (["Intercept_mean"] + [f"{d}_mean" for d in LEVEL_DIFFS]
                    + ["Intercept_std"] + [f"{d}_std" for d in LEVEL_DIFFS])
    cent = (summ[[f"{d}_mean" for d in LEVEL_DIFFS]]
            .apply(lambda x: x - x.mean(), axis=0))

    derivs = {}
    for contrast in CONTRASTS:
        pre = RELMAN / f"{contrast}_derivative_quantiles.tsv"
        if pre.exists():
            derivs[contrast] = pd.read_table(pre, sep="\t")
    if not derivs and (RELMAN / "beta_var.nc").exists():
        import xarray as xr
        beta = xr.open_dataset(RELMAN / "beta_var.nc").stack(sample=["chain", "draw"])
        names = dict(zip(beta.coords["covariate"].values, LEVEL_DIFFS))
        for contrast in CONTRASTS:
            v = cent[f"{contrast}_mean"].sort_values()
            d = (beta["beta_var"].sel(feature=v.tail(N_LOGRATIO).index).mean("feature")
                 - beta["beta_var"].sel(feature=v.head(N_LOGRATIO).index).mean("feature"))
            derivs[contrast] = (d.to_dataframe().reset_index()
                                .assign(covariate=lambda x: x["covariate"].map(names))
                                .groupby("covariate")["beta_var"]
                                .quantile([0.05, 0.5, 0.95]).reset_index()
                                .rename(columns={"level_1": "quantile"}))

    return {
        "cent": cent,
        "table": biom.load_table(str(RELMAN_TBL)).to_dataframe(dense=True).T,
        "md": pd.read_table(RELMAN_MD, sep="\t", index_col=0),
        "derivs": derivs,
    }


def style_b(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color(RULE)
        spine.set_linewidth(0.7)
    ax.tick_params(labelsize=7.5, length=3, color=RULE)
    ax.set_xlim(-0.25, len(LEVELS) - 0.75)


def panel_b_placeholder(axes):
    for ax in axes:
        ax.set_axis_off()
    ax = axes[0]
    ax.add_patch(FancyBboxPatch(
        (0.02, -1.18), 2.14, 2.16, boxstyle="round,pad=0,rounding_size=0.02",
        transform=ax.transAxes, facecolor="#F5F6F8", edgecolor=RULE,
        linewidth=0.7, linestyle=(0, (3, 2)), clip_on=False))
    ax.text(1.09, 0.20, "Antibiotic log-ratio dynamics", transform=ax.transAxes,
            ha="center", va="center", fontsize=10, fontweight="bold", color=INK)
    ax.text(1.09, -0.08,
            "awaiting BIRDMAn model outputs\n\n"
            "results/relman_abx/birdman_results.beta_var.tsv\n"
            "results/relman_abx/beta_var.nc  (or *_derivative_quantiles.tsv)",
            transform=ax.transAxes, ha="center", va="center", fontsize=7.5,
            color=MUTED, linespacing=1.7)


def panel_b(axes, data):
    (lr1, lr2), (d1, d2) = axes[:2], axes[2:]
    cent, tbl, md = data["cent"], data["table"], data["md"]
    subjects = sorted(md["host_subject_id"].unique())
    colors = dict(zip(subjects, SUBJ_COLORS))

    for ax, contrast in zip((lr1, lr2), CONTRASTS):
        v = cent[f"{contrast}_mean"].sort_values()
        v = v.loc[[f for f in v.index if f in tbl.columns]]
        num = tbl[v.tail(N_LOGRATIO).index].sum(axis=1) + 1
        den = tbl[v.head(N_LOGRATIO).index].sum(axis=1) + 1
        df = pd.DataFrame({"lr": np.log(num / den)}).join(md)
        df["x"] = df["antibiotic"].map({lvl: i for i, lvl in enumerate(LEVELS)})

        for subj, g in df.groupby("host_subject_id"):
            a = g.groupby("x")["lr"].agg(["mean", "std"]).sort_index()
            ax.plot(a.index, a["mean"], marker="o", ms=4, lw=1.3,
                    color=colors[subj], zorder=3)
            ax.fill_between(a.index, a["mean"] - a["std"], a["mean"] + a["std"],
                            color=colors[subj], alpha=0.20, lw=0, zorder=2)

        ax.set_title(contrast.replace("_vs_", " vs. "), fontsize=10, pad=6)
        ax.set_xticks(range(len(LEVELS)))
        ax.set_xticklabels([])
        style_b(ax)

    lr2.tick_params(labelleft=False)
    lr1.set_ylabel(r"$\ln\left(\frac{\sum_{\mathrm{Top\ 40\ OTUs}}}"
                   r"{\sum_{\mathrm{Bottom\ 40\ OTUs}}}\right)$", fontsize=11)
    lr1.legend(handles=[Line2D([], [], color=colors[s], lw=1.6, label=str(s))
                        for s in subjects],
               title="Subject", title_fontsize=7.5, fontsize=7, loc="upper right",
               frameon=True, framealpha=0.95, edgecolor=RULE, labelspacing=0.3)

    for ax, contrast in zip((d1, d2), CONTRASTS):
        q = data["derivs"].get(contrast)
        if q is None:
            ax.set_axis_off()
            continue
        val = "beta_var" if "beta_var" in q.columns else q.columns[-1]
        q = q.assign(covariate=pd.Categorical(q["covariate"],
                                              categories=LEVEL_DIFFS, ordered=True))
        q = q.sort_values("covariate")
        lo, mid, hi = (q[q["quantile"] == p][val].to_numpy()
                       for p in (0.05, 0.5, 0.95))
        xs = np.arange(len(mid))

        ax.fill_between(xs, lo, hi, color="#999999", alpha=0.25, lw=0, zorder=2)
        for edge in (lo, hi):
            ax.plot(xs, edge, ls="--", lw=0.9, color=MUTED, zorder=3)
        ax.plot(xs, mid, marker="o", ms=4, lw=1.3, color="#555555", zorder=4)
        ax.set_xticks(xs)
        ax.set_xticklabels(LEVELS[1:], rotation=45, ha="right", fontsize=7.5)
        style_b(ax)

    d2.tick_params(labelleft=False)
    d1.set_ylabel(r"$\Delta \ln\left(\frac{\bar{\beta}_{\mathrm{Top\ 40\ OTUs}}}"
                  r"{\bar{\beta}_{\mathrm{Bottom\ 40\ OTUs}}}\right)$", fontsize=11)


# ---------------------------------------------------------------------------
def main():
    fig = plt.figure(figsize=(13.0, 6.0))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.30, 1.0, 1.0],
                  hspace=0.12, wspace=0.10, left=0.155, right=0.985,
                  top=0.92, bottom=0.15)
    b_axes = [fig.add_subplot(gs[r, c]) for r in (0, 1) for c in (1, 2)]

    panel_a(fig.add_subplot(gs[:, 0]))

    data = load_relman()
    if data is None:
        print(f"[panel b] missing outputs under {RELMAN.relative_to(REPO)}"
              " - drawing placeholder")
        panel_b_placeholder(b_axes)
    else:
        panel_b(b_axes, data)

    for x, label in ((0.012, "a"), (0.520, "b")):
        fig.text(x, 0.945, label, fontsize=20, fontweight="bold", color=INK)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        out = OUTDIR / f"figure2.{ext}"
        fig.savefig(out, dpi=450 if ext == "png" else None, facecolor="white")
        print("wrote", out.relative_to(REPO))
    plt.close(fig)


if __name__ == "__main__":
    main()
