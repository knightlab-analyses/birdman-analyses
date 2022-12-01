#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/relman_abx/%x.out
#SBATCH --partition=short
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

import glob
import re

import arviz as az
import biom
import pandas as pd
import xarray as xr

feat_regex = re.compile("F\d{4}_(.*)\.nc")

tbl = biom.load_table("data/relman_abx/processed/processed_tbl.biom")

concatenation_name = "feature"
coords = {"feature": tbl.ids(axis="observation")}

inf_dir = "/panfs/grahman/birdman-analyses-final/relman_abx/inferences"
inf_file_list = sorted(glob.glob(f"{inf_dir}/*.nc"))

results = []
feat_names = []

for inf_file in inf_file_list:
    feat_name = feat_regex.search(inf_file).groups()[0]
    inf = az.from_netcdf(inf_file)
    post = inf.posterior
    beta_var = post["beta_var"].drop_sel(covariate=["Intercept"])
    results.append(beta_var)
    feat_names.append(feat_name)

full_beta = xr.concat(
    results,
    pd.Index(feat_names, name="feature")
)

full_beta.to_netcdf("results/relman_abx/beta_var.nc")
