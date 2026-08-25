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

Panel a runs from data committed to this repository. Panel b additionally needs
the BIRDMAn model outputs under results/relman_abx/, which are not committed;
when they are absent a labelled placeholder is drawn so the composite still
renders.

HDI unpacking follows https://github.com/ahdilmore/MARS_Birdman
(birdman/zebra_birdman_analysis.ipynb): split on the comma, strip the
parentheses, flag credible when the interval excludes zero, then convert the
bounds to offsets from the mean for errorbar().

Usage:
    python scripts/F2_figure2.py [-n 10] [--include-non-credible] [--no-ci]
"""

import argparse
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
N_LOGRATIO = 40

# Panel a: two sequential ramps about a pale midpoint (a diverging scheme keyed
# to magnitude). Panel b subjects: categorical, validated colourblind-safe --
# chroma >= 0.1, adjacent CVD dE 10.5 (protan), normal-vision dE 26.9.
BV_CMAP, HEALTH_CMAP = plt.get_cmap("Reds"), plt.get_cmap("Blues")
RAMP_LO, RAMP_HI = 0.22, 0.92
BV_KEY, HEALTH_KEY = BV_CMAP(0.75), HEALTH_CMAP(0.55)
SUBJ_COLORS = ["#0077BB", "#EE7733", "#009988"]
INK, MUTED, RULE, GRID = "#222222", "#666666", "#BBBBBB", "#DDDDDD"


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
                    fmt="none", ecolor=INK, elinewidth=0.6, capsize=1.5,
                    capthick=0.6, alpha=0.8, zorder=3)
    ax.axvline(0, color=INK, linewidth=0.8, zorder=4)

    ax.set_yticks(ys)
    ax.set_yticklabels([f.replace("_", " ") for f in sel.index], fontsize=8)
    ax.set_ylim(-0.7, len(sel) - 0.3)
    ax.set_xlabel("BIRDMAn Differential", fontsize=9.5, color=INK, labelpad=6)
    ax.tick_params(axis="both", labelsize=8, length=3, color=RULE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(RULE)

    # Legend below the panel, horizontal — as in the reference layout.
    ax.legend(handles=[Patch(facecolor=BV_KEY, label="BV-associated"),
                       Patch(facecolor=HEALTH_KEY, label="Healthy-associated")],
              loc="upper center", bbox_to_anchor=(0.5, -0.10), ncol=2,
              fontsize=9, frameon=False, handlelength=1.4, handleheight=1.0,
              columnspacing=2.0, borderpad=0.0)


# ---------------------------------------------------------------------------
# Panel b — antibiotic log-ratio dynamics
# ---------------------------------------------------------------------------
def log_ratio(table, top_feats, bot_feats):
    lr = pd.concat([table.loc[:, top_feats].sum(axis=1),
                    table.loc[:, bot_feats].sum(axis=1)], axis=1)
    lr.columns = ["num", "denom"]
    lr = lr.dropna(how="all") + 1
    return np.log(lr["num"] / lr["denom"])


def load_relman():
    """Return the pieces panel b needs, or None if the model outputs are absent."""
    beta_tsv = RELMAN / "birdman_results.beta_var.tsv"
    if not (beta_tsv.exists() and RELMAN_TBL.exists() and RELMAN_MD.exists()):
        return None

    import biom

    summ = pd.read_table(beta_tsv, sep="\t", index_col=0)
    summ.index = summ.index.astype(str)
    # Columns arrive as Intercept, 6 means, Intercept_std, 6 stds — positional,
    # matching scripts/relman_abx/2.01-summarize_birdman.py.
    rename = dict(zip(summ.columns[1:7], [f"{d}_mean" for d in LEVEL_DIFFS]))
    rename.update(dict(zip(summ.columns[8:], [f"{d}_std" for d in LEVEL_DIFFS])))
    summ = summ.rename(columns=rename)

    cent = (summ[[f"{d}_mean" for d in LEVEL_DIFFS]]
            .apply(lambda x: x - x.mean(), axis=0))
    tbl = biom.load_table(str(RELMAN_TBL)).to_dataframe(dense=True).T
    md = pd.read_table(RELMAN_MD, sep="\t", index_col=0)

    derivs = {}
    for contrast in CONTRASTS:
        pre = RELMAN / f"{contrast}_derivative_quantiles.tsv"
        if pre.exists():
            derivs[contrast] = pd.read_table(pre, sep="\t")
    if not derivs and (RELMAN / "beta_var.nc").exists():
        import xarray as xr
        beta = xr.open_dataset(RELMAN / "beta_var.nc").stack(sample=["chain", "draw"])
        cmap = dict(zip(beta.coords["covariate"].values, LEVEL_DIFFS))
        for contrast in CONTRASTS:
            vals = cent[f"{contrast}_mean"].sort_values()
            diff = (beta["beta_var"].sel(feature=vals.tail(N_LOGRATIO).index).mean("feature")
                    - beta["beta_var"].sel(feature=vals.head(N_LOGRATIO).index).mean("feature"))
            derivs[contrast] = (diff.to_dataframe().reset_index()
                                .assign(covariate=lambda x: x["covariate"].map(cmap))
                                .groupby("covariate")["beta_var"]
                                .quantile([0.05, 0.5, 0.95]).reset_index()
                                .rename(columns={"level_1": "quantile"}))
    return {"cent": cent, "table": tbl, "md": md, "derivs": derivs}


def style_b_axis(ax):
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color(RULE)
        s.set_linewidth(0.7)
    ax.tick_params(labelsize=7.5, length=3, color=RULE)


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


def panel_b(axes, data):
    (lr1, lr2), (d1, d2) = axes[:2], axes[2:]
    cent, tbl, md = data["cent"], data["table"], data["md"]
    order = {lvl: i for i, lvl in enumerate(LEVELS)}
    subjects = sorted(md["host_subject_id"].unique())
    colors = dict(zip(subjects, SUBJ_COLORS))

    for ax, contrast in zip((lr1, lr2), CONTRASTS):
        vals = cent[f"{contrast}_mean"].sort_values()
        vals = vals.loc[[f for f in vals.index if f in tbl.columns]]
        lr = log_ratio(tbl, vals.tail(N_LOGRATIO).index, vals.head(N_LOGRATIO).index)
        df = pd.DataFrame({"log_ratio": lr}).join(md[["antibiotic", "host_subject_id"]])
        df["x"] = df["antibiotic"].map(order)
        for subj, g in df.groupby("host_subject_id"):
            agg = g.groupby("x")["log_ratio"].agg(["mean", "std"]).sort_index()
            ax.plot(agg.index, agg["mean"], marker="o", ms=4, lw=1.3,
                    color=colors[subj], zorder=3)
            ax.fill_between(agg.index, agg["mean"] - agg["std"],
                            agg["mean"] + agg["std"], color=colors[subj],
                            alpha=0.20, linewidth=0, zorder=2)
        ax.set_title(contrast.replace("_vs_", " vs. "), fontsize=10, pad=6)
        ax.set_xticks(range(len(LEVELS)))
        ax.set_xticklabels([])
        ax.set_xlim(-0.25, len(LEVELS) - 0.75)
        style_b_axis(ax)

    lr2.tick_params(labelleft=False)
    lr1.set_ylabel(r"$\ln\left(\frac{\sum_{\mathrm{Top\ 40\ OTUs}}}"
                   r"{\sum_{\mathrm{Bottom\ 40\ OTUs}}}\right)$", fontsize=11)
    lr1.legend(handles=[Line2D([], [], color=colors[s], lw=1.6, label=str(s))
                        for s in subjects],
               title="Subject", title_fontsize=7.5, fontsize=7,
               loc="upper right", frameon=True, framealpha=0.95,
               edgecolor=RULE, labelspacing=0.3, handlelength=1.6)

    for ax, contrast in zip((d1, d2), CONTRASTS):
        q = data["derivs"].get(contrast)
        if q is None:
            ax.set_axis_off()
            continue
        q = q.copy()
        q["covariate"] = pd.Categorical(q["covariate"], categories=LEVEL_DIFFS,
                                        ordered=True)
        q = q.sort_values("covariate")
        val = "beta_var" if "beta_var" in q.columns else q.columns[-1]
        lo = q[q["quantile"] == 0.05][val].to_numpy()
        mid = q[q["quantile"] == 0.50][val].to_numpy()
        hi = q[q["quantile"] == 0.95][val].to_numpy()
        xs = np.arange(len(mid))
        ax.fill_between(xs, lo, hi, color="#999999", alpha=0.25, linewidth=0,
                        zorder=2)
        for edge in (lo, hi):
            ax.plot(xs, edge, ls="--", lw=0.9, color=MUTED, zorder=3)
        ax.plot(xs, mid, marker="o", ms=4, lw=1.3, color="#555555", zorder=4)
        ax.set_xticks(xs)
        ax.set_xticklabels(LEVELS[1:], rotation=45, ha="right", fontsize=7.5)
        ax.set_xlim(-0.25, len(xs) - 0.75)
        style_b_axis(ax)

    d2.tick_params(labelleft=False)
    d1.set_ylabel(r"$\Delta \ln\left(\frac{\bar{\beta}_{\mathrm{Top\ 40\ OTUs}}}"
                  r"{\bar{\beta}_{\mathrm{Bottom\ 40\ OTUs}}}\right)$", fontsize=11)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=10, help="features per side in panel a")
    ap.add_argument("--include-non-credible", action="store_true",
                    help="do not require the HDI to exclude zero in panel a")
    ap.add_argument("--no-ci", action="store_true",
                    help="omit HDI error bars from panel a (reference styling)")
    args = ap.parse_args()

    fig = plt.figure(figsize=(13.0, 6.0))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.30, 1.0, 1.0],
                  hspace=0.12, wspace=0.10, left=0.155, right=0.985,
                  top=0.92, bottom=0.15)
    ax_a = fig.add_subplot(gs[:, 0])
    b_axes = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]),
              fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])]

    panel_a(ax_a, args.n, credible_only=not args.include_non_credible,
            show_ci=not args.no_ci)

    data = load_relman()
    if data is None:
        print(f"[panel b] missing model outputs under {RELMAN.relative_to(REPO)}"
              " - drawing placeholder")
        panel_b_placeholder(b_axes)
    else:
        panel_b(b_axes, data)

    # Panel labels in figure coordinates so long tick labels cannot push them
    # off-canvas.
    for x, label in ((0.012, "a"), (0.520, "b")):
        fig.text(x, 0.945, label, fontsize=20, fontweight="bold",
                 va="bottom", ha="left", color=INK)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        out = OUTDIR / f"figure2.{ext}"
        fig.savefig(out, dpi=450 if ext == "png" else None, facecolor="white")
        print("wrote", out.relative_to(REPO))
    plt.close(fig)


if __name__ == "__main__":
    main()
