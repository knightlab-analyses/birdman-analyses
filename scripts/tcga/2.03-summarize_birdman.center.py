#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

import glob
import re

import arviz as az
import pandas as pd

FEAT_REGEX = re.compile("F\d{4}_(.*).nc")

inference_dir = "/panfs/grahman/birdman-analyses-final/tcga/species/bacteria"
outfile = "results/tcga/species/birdman_results.bacteria.center.nc"
all_inf_files = glob.glob(f"{inference_dir}/*.nc")

feat_inf_list = []
feat_names = []
for inf_file in all_inf_files:
    this_feat_id = FEAT_REGEX.search(inf_file).groups()[0]
    this_feat_inf = az.from_netcdf(inf_file).posterior["center_re"]

    feat_inf_list.append(this_feat_inf)
    feat_names.append(this_feat_id)

full_inf = xr.concat(
    feat_inf_list,
    pd.Index(feat_names, name="feature"),
)
full_inf.to_netcdf(outfile)
