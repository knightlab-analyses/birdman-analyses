#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/sample_size/%x.out
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

import glob
import os
import re

import arviz as az
import pandas as pd
import xarray as xr

inf_files = glob.glob("results/sample_size/sim/*.nc")

post_dict = dict()
for f in inf_files:
    sample_size = int(re.search("(\d+).nc", os.path.basename(f)).groups()[0])
    post = az.from_netcdf(f).posterior
    post_dict[sample_size] = az.from_netcdf(f).posterior

order = list(map(str, sorted(post_dict.keys())))

beta_var_dfs = []
inv_disp_dfs = []
for sample_size, post in post_dict.items():
    _df = (
        post["beta_var"]
        .stack(sample=["chain", "draw"])
        .to_dataframe()
        .reset_index()
        .assign(sample_size=str(sample_size))
        .melt(id_vars=["sample_size", "column"], value_vars=["beta_var"])
        .drop(columns=["variable"])
    )
    beta_var_dfs.append(_df)

    _df = (
        post["inv_disp"]
        .stack(sample=["chain", "draw"])
        .to_dataframe()
        .reset_index()
        .assign(sample_size=str(sample_size))
        .melt(id_vars=["sample_size"], value_vars=["inv_disp"])
        .drop(columns=["variable"])
    )
    inv_disp_dfs.append(_df)

beta_var_df = pd.concat(beta_var_dfs).reset_index(drop=True)
beta_var_df["sample_size"] = pd.Categorical(beta_var_df["sample_size"],
                                            categories=order,
                                            ordered=True)
beta_var_df.to_csv("results/sample_size/post_beta_var_df.tsv", sep="\t")

inv_disp_df = pd.concat(inv_disp_dfs).reset_index(drop=True)
inv_disp_df["sample_size"] = pd.Categorical(inv_disp_df["sample_size"],
                                            categories=order,
                                            ordered=True)
inv_disp_df.to_csv("results/sample_size/post_inv_disp_df.tsv", sep="\t")
