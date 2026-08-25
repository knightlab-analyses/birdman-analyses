#!/usr/bin/env python
"""collapse the concatenated posterior into the two small tables that
notebooks/relman_abx_figures.ipynb needs for the bottom row of figure 2b.

beta_var.nc is ~79 MB, too big to commit, but the figure only ever consumes
6 covariates x 3 quantiles from it. run this once after 2.03-concat_infs.py
and commit the output instead.
"""
import os

import pandas as pd
import xarray as xr

PROJ_DIR = os.environ.get(
    "BIRDMAN_ANALYSES_DIR",
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)
RESULTS_DIR = f"{PROJ_DIR}/results/relman_abx"

levels = ["preCp", "FirstCp", "FirstWPC", "Interim", "SecondCp", "SecondWPC",
          "PostCp"]
levels_diffs = [levels[i] + "_vs_" + levels[i - 1]
                for i in range(1, len(levels))]

# same positional rename the notebook does: patsy names the backward-difference
# columns after the *first* level of each pair, so order is load-bearing
summ_diff = pd.read_table(f"{RESULTS_DIR}/birdman_results.beta_var.tsv",
                          sep="\t", index_col=0)
summ_diff.index = summ_diff.index.astype(str)
rename_dict = dict(zip(summ_diff.columns[1:7],
                       [x + "_mean" for x in levels_diffs]))
rename_dict.update(dict(zip(summ_diff.columns[8:],
                            [x + "_std" for x in levels_diffs])))
summ_diff = summ_diff.rename(columns=rename_dict)

mean_level_diffs = summ_diff.columns[1:7]
summ_diff_cent = summ_diff[mean_level_diffs].apply(lambda x: x - x.mean(),
                                                   axis=0)

beta_var = xr.open_dataset(f"{RESULTS_DIR}/beta_var.nc")
abx_contrast_dict = dict(zip(beta_var.coords["covariate"].values, levels_diffs))


def get_derivative(vals, n=40):
    _vals = vals.sort_values()
    _top = _vals.tail(n).index
    _bot = _vals.head(n).index
    per_draw_top = beta_var.sel({"feature": _top}).mean(["feature"])
    per_draw_bot = beta_var.sel({"feature": _bot}).mean(["feature"])

    return (
        (per_draw_top - per_draw_bot)
        .to_dataframe()
        .reset_index()
        .assign(covariate=lambda x: x["covariate"].map(abx_contrast_dict))
        .groupby("covariate")["beta_var"]
        .quantile([0.05, 0.5, 0.95])
        .rename_axis(index=["covariate", "quantile"])
        .reset_index()
    )


for col, prefix in [("FirstCp_vs_preCp_mean", "fcp"),
                    ("SecondCp_vs_Interim_mean", "scp")]:
    out = get_derivative(summ_diff_cent[col])
    outfile = f"{RESULTS_DIR}/{prefix}_derivative_quantiles.tsv"
    out.to_csv(outfile, sep="\t", index=False)
    print(f"wrote {outfile}")
