#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/sampling_depth/%x.out
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

inf_files = glob.glob("results/sampling_depth/sim/*.nc")

post_dict = dict()
for f in inf_files:
    mean_depth = int(re.search("(\d+).nc", os.path.basename(f)).groups()[0])
    post = az.from_netcdf(f).posterior
    post_dict[mean_depth] = az.from_netcdf(f).posterior

concat_post = xr.concat(
    [post for post in post_dict.values()],
    dim=pd.Index(post_dict.keys(), name="mean_depth")
)
concat_post.to_netcdf("results/sampling_depth/concat_posterior.nc")
