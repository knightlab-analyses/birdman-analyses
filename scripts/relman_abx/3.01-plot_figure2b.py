#!/usr/bin/env python
"""figure 2b -- antibiotic-course dynamics for the relman_abx case study.

reads only committed artifacts, so this runs from a fresh clone with no cluster
and no model refit:
    results/relman_abx/birdman_results.beta_var.tsv
    results/relman_abx/{fcp,scp}_derivative_quantiles.tsv
    data/relman_abx/processed/{processed_tbl.biom,processed_md.tsv}

standalone:
    python scripts/relman_abx/3.01-plot_figure2b.py [-o out.pdf]

embedded in a combined figure (e.g. alongside 2a):
    from importlib.machinery import SourceFileLoader
    fig2b = SourceFileLoader("fig2b", ".../3.01-plot_figure2b.py").load_module()

    fig = plt.figure(figsize=(11, 14))
    outer = GridSpec(nrows=2, ncols=1, height_ratios=[1, 2])
    my_panel_a(fig, outer[0])
    fig2b.plot_figure_2b(fig=fig, subplot_spec=outer[1])

note: the top row's shaded bands come from seaborn's bootstrap CI, whose RNG is
unseeded by default -- pass seed=<int> for run-to-run identical output.
"""
import argparse
import os

import biom
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Patch

PROJ_DIR = os.environ.get(
    "BIRDMAN_ANALYSES_DIR",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

LEVELS = ["preCp", "FirstCp", "FirstWPC", "Interim", "SecondCp", "SecondWPC",
          "PostCp"]
LEVELS_DIFFS = [LEVELS[i] + "_vs_" + LEVELS[i - 1] for i in range(1, len(LEVELS))]
SUBJECTS = ["494.D", "494.E", "494.F"]
# the two contrasts the figure selects features on
CONTRASTS = ["FirstCp_vs_preCp", "SecondCp_vs_Interim"]

TOP_LABEL = (r"$\ln\left(\frac{\sum\ \mathrm{Top\ 40\ OTUs}}"
             r"{\sum \mathrm{Bottom\ 40\ OTUs}}\right)$")
BOT_LABEL = (r"$\Delta\ \ln\left(\frac {\bar{\beta}_{\mathrm{Top\ 40\ OTUs}}}"
             r" {\bar{\beta}_{\mathrm{Bottom\ 40\ OTUs}}}\right)$")


def _log_ratio(table, md, top_feats, bot_feats):
    lr = pd.concat([table.loc[:, top_feats].sum(axis=1),
                    table.loc[:, bot_feats].sum(axis=1)], axis=1)
    lr.columns = ["num", "denom"]
    lr = lr.dropna(how="all") + 1
    lr["log_ratio"] = np.log(lr["num"] / lr["denom"])
    return lr.join(md)


def load_data(proj_dir=PROJ_DIR, n=40):
    """committed artifacts -> everything the panels need."""
    res = f"{proj_dir}/results/relman_abx"
    dat = f"{proj_dir}/data/relman_abx/processed"

    summ = pd.read_table(f"{res}/birdman_results.beta_var.tsv", index_col=0)
    summ.index = summ.index.astype(str)
    # patsy names backward-difference columns after the first level of each
    # pair, so the rename is positional -- column order is load-bearing
    ren = dict(zip(summ.columns[1:7], [x + "_mean" for x in LEVELS_DIFFS]))
    ren.update(dict(zip(summ.columns[8:], [x + "_std" for x in LEVELS_DIFFS])))
    summ = summ.rename(columns=ren)
    mean_cols = summ.columns[1:7]
    summ_cent = summ[mean_cols].apply(lambda x: x - x.mean(), axis=0)

    table = biom.load_table(f"{dat}/processed_tbl.biom").to_dataframe(dense=True).T
    md = pd.read_table(f"{dat}/processed_md.tsv", index_col=0)
    md["antibiotic"] = pd.Categorical(md["antibiotic"], categories=LEVELS,
                                      ordered=True)

    predictors = dict()
    for diff in LEVELS_DIFFS:
        ranked = summ_cent[f"{diff}_mean"].sort_values(ascending=False)
        predictors[diff] = _log_ratio(
            table, md, ranked.head(n).index, ranked.tail(n).index
        )["log_ratio"]
    predictor_df = pd.DataFrame.from_dict(predictors).join(md)

    derivs = {}
    for key in CONTRASTS:
        d = pd.read_table(f"{res}/{key}_derivative_quantiles.tsv")
        d["covariate"] = pd.Categorical(d["covariate"], categories=LEVELS_DIFFS,
                                        ordered=True)
        derivs[key] = d.sort_values(by="covariate")

    return {"predictor_df": predictor_df, "derivs": derivs}


def plot_figure_2b(fig=None, subplot_spec=None, data=None, seed=None):
    """draw figure 2b; returns (fig, axes).

    fig/subplot_spec let a caller embed this inside a larger combined figure.
    """
    if data is None:
        data = load_data()

    if fig is None:
        fig = plt.figure(figsize=(11, 10))
    if subplot_spec is None:
        gs = GridSpec(ncols=2, nrows=3, figure=fig)
    else:
        # the two columns share a y-axis and hide the right labels, so keep them
        # tight -- the caller controls the gap to whatever sits left of them
        gs = GridSpecFromSubplotSpec(ncols=2, nrows=2, subplot_spec=subplot_spec,
                                     wspace=0.08, hspace=0.18)

    palette = dict(zip(SUBJECTS, sns.color_palette("colorblind", 3)))
    patches = [Patch(facecolor=c, label=s, linewidth=0.5, edgecolor="black")
               for s, c in palette.items()]

    lr1 = fig.add_subplot(gs[0, 0])
    lr2 = fig.add_subplot(gs[0, 1], sharey=lr1)
    args = {"data": data["predictor_df"], "x": "antibiotic",
            "hue": "host_subject_id", "marker": "o", "palette": palette,
            "legend": False}
    if seed is not None:
        args["seed"] = seed
    for ax, contrast in zip((lr1, lr2), CONTRASTS):
        sns.lineplot(**args, y=contrast, ax=ax)
        ax.set_xticklabels([])
        ax.set_xlabel("")
        ax.grid()
        ax.set_xlim([-0.25, 6.25])
    lr1.set_ylabel(TOP_LABEL, fontsize="x-large")
    lr2.set_ylabel("")
    lr2.tick_params("y", width=0, labelleft=False)
    lr1.set_title("FirstCp vs. preCp")
    lr2.set_title("SecondCp vs. Interim")
    lr1.legend(handles=patches, title="Subject")

    d1 = fig.add_subplot(gs[1, 0])
    d2 = fig.add_subplot(gs[1, 1], sharey=d1)
    for ax, contrast in zip((d1, d2), CONTRASTS):
        _df = data["derivs"][contrast]
        lo = _df.query("quantile == 0.05")["beta_var"]
        mid = _df.query("quantile == 0.5")["beta_var"]
        hi = _df.query("quantile == 0.95")["beta_var"]
        x = np.arange(6)
        ax.fill_between(x, lo, hi, color="gray", alpha=0.25)
        ax.plot(x, mid, markeredgecolor="black", markeredgewidth=0.5,
                color="gray", marker="o")
        for bound in (hi, lo):
            ax.plot(x, bound, color="gray", zorder=1, ls="--")
        ax.grid()
        ax.set_xticks(np.arange(-1, 6))
        ax.set_xticklabels(LEVELS, rotation=45, ha="right", fontsize="large")
        ax.set_xlim([-1.25, 5.25])
    d1.set_ylabel(BOT_LABEL, fontsize="x-large")
    d2.tick_params("y", width=0, labelleft=False)

    return fig, {"lr": (lr1, lr2), "deriv": (d1, d2)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-o", "--outfile",
                   default=f"{PROJ_DIR}/figures/relman_abx/abx_dynamics.pdf")
    p.add_argument("--seed", type=int, default=None,
                   help="seed seaborn's bootstrap CI for reproducible bands")
    args = p.parse_args()

    plt.style.use(f"{PROJ_DIR}/notebooks/paper.mplstyle")
    fig, _ = plot_figure_2b(seed=args.seed)
    os.makedirs(os.path.dirname(args.outfile), exist_ok=True)
    fig.savefig(args.outfile)
    print(f"wrote {args.outfile}")


if __name__ == "__main__":
    main()
