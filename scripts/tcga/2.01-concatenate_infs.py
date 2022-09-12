#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga/%x.out
#SBATCH --partition=short
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

import glob
import json

import arviz as az
import biom
from birdman.model_util import concatenate_inferences
import numpy as np
import xarray as xr

tbl = biom.load_table("data/tcga/processed/processed_tbl.fungi.biom")

concatenation_name = "feature"
coords = {"feature": tbl.ids(axis="observation")}

inf_dir = "/panfs/grahman/birdman-analyses-final/tcga/WIS_bacteria_fungi"
inf_file_list = sorted(glob.glob(f"{inf_dir}/*.nc"))
inf_list = [az.from_netcdf(x) for x in inf_file_list]

all_inf = concatenate_inferences(
    inf_list,
    coords=coords,
    concatenation_name="feature"
)
all_inf.to_netcdf("results/tcga/inf.fungi.nc")
