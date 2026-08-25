"""
Figure 2 — BIRDMAn use cases, rendered in full.

    a  Vaginal microbiome (Ravel et al. 2011): BIRDMAn posterior differentials
       with 95% HDI error bars, for the N most BV-associated and N most
       health-associated features. Only features whose HDI excludes zero are
       eligible; ranking on the posterior mean alone promotes features whose
       means are large only because they are barely estimated.

    b  Dual-course ciprofloxacin (Dethlefsen & Relman 2011): per-subject sample
       log-ratio trajectories for the FirstCp and SecondCp contrasts (top row),
       and the corresponding posterior derivative with a 5-95% ribbon (bottom).

Panel a runs from data committed to this repository. Panel b additionally needs
the BIRDMAn model outputs under results/relman_abx/, which are not committed
(see repo notes); when they are absent a labelled placeholder is drawn in their
place so the composite still renders.

HDI unpacking follows https://github.com/ahdilmore/MARS_Birdman
(birdman/zebra_birdman_analysis.ipynb): split on the comma, strip the
parentheses, flag credible when the interval excludes zero, then convert the
bounds to offsets from the mean for errorbar().

Usage:
    python scripts/F2_figure2.py [-n 10] [--include-non-credible]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

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

# Diverging pair (polarity) and categorical subjects — both validated
# colourblind-safe: chroma >= 0.1, adjacent CVD dE >= 8, normal-vision dE >= 15.
BV_COLOR, HEALTH_COLOR = "#CC3311", "#0077BB"
SUBJ_COLORS = ["#0077BB", "#EE7733", "#009988"]
INK, MUTED, RULE = "#222222", "#666666", "#BBBBBB"


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


def italicize(name):
    if name.startswith("Lactobacillus_"):
        return f"L. {name.split('_', 1)[1]}", True
    label = name.replace("_", " ")
    return label, len(label.split()) > 1


def panel_a(ax, n, credible_only):
    df = (pd.read_csv(DIFFERENTIALS, sep="\t")
            .rename(columns={"feature id": "Feature"})
            .set_index("Feature"))
    df = unpack_hdi(df, COVARIATE)
    mean_col = f"{COVARIATE}_mean"
    n_total = len(df)
    if credible_only:
        df = df[df["credible"] == "yes"]

    df = df.sort_values(mean_col)
    sel = df if len(df) <= 2 * n else pd.concat([df.head(n), df.tail(n)])

    ys = np.arange(len(sel))
    means = sel[mean_col].to_numpy()
    colors = [BV_COLOR if m > 0 else HEALTH_COLOR for m in means]

    ax.barh(ys, means, color=colors, height=0.68, zorder=2,
            edgecolor="white", linewidth=0.4)
    ax.errorbar(means, ys, xerr=sel[["lower", "upper"]].to_numpy().T, fmt="none",
                ecolor=INK, elinewidth=0.7, capsize=1.8, capthick=0.7, zorder=3)
    ax.axvline(0, color=RULE, linewidth=0.8, zorder=1)

    labels, italics = zip(*(italicize(f) for f in sel.index))
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=6.5)
    for tick, is_it in zip(ax.get_yticklabels(), italics):
        if is_it:
            tick.set_style("italic")
    ax.set_ylim(-0.7, len(sel) - 0.3)
    ax.set_xlabel("BIRDMAn differential\n(log fold change, BV vs healthy)",
                  fontsize=7.5, color=INK)
    ax.tick_params(axis="x", labelsize=6.5)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(RULE)

    ax.legend(handles=[Line2D([], [], color=BV_COLOR, lw=4, label="BV-associated"),
                       Line2D([], [], color=HEALTH_COLOR, lw=4, label="Health-associated")],
              loc="upper left", fontsize=6, frameon=False,
              handlelength=1.0, borderpad=0.2, labelspacing=0.3)

    note = f"{len(sel)} of {n_total} features · 95% HDI"
    if credible_only:
        note += " · HDI excludes 0"
    ax.text(0.0, 1.005, note, transform=ax.transAxes, fontsize=6,
            color=MUTED, va="bottom")


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
    # Columns arrive as Intercept, 6 means, Intercept_std, 6 stds - positional,
    # matching scripts/relman_abx/2.01-summarize_birdman.py.
    rename = dict(zip(summ.columns[1:7], [f"{d}_mean" for d in LEVEL_DIFFS]))
    rename.update(dict(zip(summ.columns[8:], [f"{d}_std" for d in LEVEL_DIFFS])))
    summ = summ.rename(columns=rename)

    mean_cols = [f"{d}_mean" for d in LEVEL_DIFFS]
    cent = summ[mean_cols].apply(lambda x: x - x.mean(), axis=0)

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
            top, bot = vals.tail(N_LOGRATIO).index, vals.head(N_LOGRATIO).index
            diff = (beta["beta_var"].sel(feature=top).mean("feature")
                    - beta["beta_var"].sel(feature=bot).mean("feature"))
            q = (diff.to_dataframe().reset_index()
                 .assign(covariate=lambda x: x["covariate"].map(cmap))
                 .groupby("covariate")["beta_var"]
                 .quantile([0.05, 0.5, 0.95]).reset_index()
                 .rename(columns={"level_1": "quantile"}))
            derivs[contrast] = q

    return {"cent": cent, "table": tbl, "md": md, "derivs": derivs}


def panel_b_placeholder(axes):
    for ax in axes:
        ax.set_axis_off()
    ax = axes[0]
    ax.add_patch(FancyBboxPatch(
        (0.02, 0.02), 1.96, 0.96, boxstyle="round,pad=0,rounding_size=0.02",
        transform=ax.transAxes, facecolor="#F2F4F7", edgecolor=RULE,
        linewidth=0.6, linestyle=(0, (3, 2)), clip_on=False, zorder=1))
    ax.text(1.0, 0.62, "Antibiotic log-ratio dynamics", transform=ax.transAxes,
            ha="center", va="top", fontsize=7.5, fontweight="bold", color=INK)
    ax.text(1.0, 0.46,
            "awaiting BIRDMAn model outputs\n\n"
            "results/relman_abx/birdman_results.beta_var.tsv\n"
            "results/relman_abx/beta_var.nc  (or *_derivative_quantiles.tsv)",
            transform=ax.transAxes, ha="center", va="top", fontsize=6,
            color=MUTED, linespacing=1.6)


def panel_b(axes, data):
    (lr1, lr2), (d1, d2) = axes[:2], axes[2:]
    cent, tbl, md = data["cent"], data["table"], data["md"]

    order = {lvl: i for i, lvl in enumerate(LEVELS)}
    subjects = sorted(md["host_subject_id"].unique())
    colors = dict(zip(subjects, SUBJ_COLORS))

    for ax, contrast in zip((lr1, lr2), CONTRASTS):
        vals = cent[f"{contrast}_mean"].sort_values()
        feats = [f for f in vals.index if f in tbl.columns]
        vals = vals.loc[feats]
        lr = log_ratio(tbl, vals.tail(N_LOGRATIO).index, vals.head(N_LOGRATIO).index)
        df = pd.DataFrame({"log_ratio": lr}).join(md[["antibiotic", "host_subject_id"]])
        df["x"] = df["antibiotic"].map(order)
        for subj, g in df.groupby("host_subject_id"):
            g = g.groupby("x", as_index=False)["log_ratio"].agg(["mean", "std"])
            ax.plot(g.index, g["mean"], marker="o", ms=2.5, lw=0.9,
                    color=colors[subj], zorder=3)
            ax.fill_between(g.index, g["mean"] - g["std"], g["mean"] + g["std"],
                            color=colors[subj], alpha=0.18, linewidth=0, zorder=2)
        ax.set_title(contrast.replace("_vs_", " vs. "), fontsize=7)
        ax.set_xticks(range(len(LEVELS)))
        ax.set_xticklabels([])
        ax.axhline(0, color=RULE, lw=0.6, ls="--", zorder=1)
        ax.tick_params(labelsize=6)

    lr1.set_ylabel(r"$\ln\left(\frac{\sum \mathrm{Top\ 40}}{\sum \mathrm{Bottom\ 40}}\right)$",
                   fontsize=7)
    lr1.legend(handles=[Line2D([], [], color=colors[s], lw=1.2, marker="o",
                               ms=2.5, label=str(s)) for s in subjects],
               title="Subject", title_fontsize=6, fontsize=6, frameon=False,
               loc="best", labelspacing=0.25, handlelength=1.4)

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
        ax.fill_between(xs, lo, hi, color=MUTED, alpha=0.20, linewidth=0, zorder=2)
        ax.plot(xs, mid, marker="o", ms=2.5, lw=0.9, color=INK, zorder=3)
        ax.axhline(0, color=RULE, lw=0.6, ls="--", zorder=1)
        ax.set_xticks(xs)
        ax.set_xticklabels([d.split("_vs_")[0] for d in LEVEL_DIFFS],
                           rotation=45, ha="right", fontsize=5.5)
        ax.tick_params(labelsize=6)
    d1.set_ylabel(r"$\Delta \ln\left(\frac{\beta_{\mathrm{Top\ 40}}}"
                  r"{\beta_{\mathrm{Bottom\ 40}}}\right)$", fontsize=7)


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=10, help="features per side in panel a")
    ap.add_argument("--include-non-credible", action="store_true",
                    help="do not require the HDI to exclude zero in panel a")
    args = ap.parse_args()

    fig = plt.figure(figsize=(9.0, 5.2))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1.15, 1.0, 1.0],
                  hspace=0.28, wspace=0.28)
    ax_a = fig.add_subplot(gs[:, 0])
    b_axes = [fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2]),
              fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[1, 2])]

    panel_a(ax_a, args.n, credible_only=not args.include_non_credible)

    data = load_relman()
    if data is None:
        print(f"[panel b] missing model outputs under {RELMAN.relative_to(REPO)}"
              " - drawing placeholder")
        panel_b_placeholder(b_axes)
    else:
        panel_b(b_axes, data)

    for ax, label in ((ax_a, "a"), (b_axes[0], "b")):
        ax.text(-0.18, 1.03, label, transform=ax.transAxes, fontsize=10,
                fontweight="bold", va="bottom", ha="right")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        out = OUTDIR / f"figure2.{ext}"
        fig.savefig(out, dpi=450 if ext == "png" else None,
                    bbox_inches="tight", facecolor="white")
        print("wrote", out.relative_to(REPO))
    plt.close(fig)


if __name__ == "__main__":
    main()
