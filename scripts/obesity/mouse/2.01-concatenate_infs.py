#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/mouse/%x.out
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

tbl = biom.load_table("data/obesity/processed/mouse/tbl_merged.genus.biom")

concatenation_name = "feature"
coords = {"feature": tbl.ids(axis="observation")}

inf_dir = "/panfs/grahman/birdman-analyses-final/obesity/mouse/inferences_genus"
inf_file_list = sorted(glob.glob(f"{inf_dir}/*.nc"))
inf_list = [az.from_netcdf(x) for x in inf_file_list]

all_inf = concatenate_inferences(
    inf_list,
    coords=coords,
    concatenation_name="feature"
)
all_inf.to_netcdf("results/obesity/mouse/inf.genus.mouse.nc")

# Save parameters for use in data-driven simulation

base_phis = (
    all_inf.posterior["base_phi"]
    .stack(sample=["chain", "draw"])
    .median("sample")
    .values
)
np.save("results/obesity/mouse/base_phis.npy", base_phis)

batch_disps = (
    all_inf.posterior["study_disp"]
    .stack(sample=["chain", "draw"])
    .median("sample")
    .values
    .ravel()
)
np.save("results/obesity/mouse/batch_disps.npy", batch_disps)

batch_offsets = (
    all_inf.posterior["study_re"]
    .stack(sample=["chain", "draw"])
    .median("sample")
    .values
    .ravel()
)
np.save("results/obesity/mouse/batch_offsets.npy", batch_offsets)
