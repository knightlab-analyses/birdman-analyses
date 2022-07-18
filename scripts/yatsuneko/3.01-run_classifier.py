#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/yatsuneko/%x.out
#SBATCH --partition=short
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

import arviz as az
import biom
import cmdstanpy
import numpy as np
import pandas as pd
from patsy import dmatrix

test_tbl = biom.load_table("data/yatsuneko/processed/processed_tbl.test.biom")
D, N = test_tbl.shape

test_md = pd.read_table("data/yatsuneko/processed/processed_md.test.tsv",
                        sep="\t", index_col=0)
dmat = pd.DataFrame(np.ones((N, 2)))

inf_train = az.from_netcdf("results/yatsuneko/inf.train.nc")
posterior = inf_train.posterior
posterior = posterior.stack(sample=["chain", "draw"])

dims = ["sample", "covariate", "feature"]
post_beta_var = posterior["beta_var"].transpose(*dims).values
post_phi = posterior["phi"].transpose("sample", "feature").values

test_model_path = "src/yatsuneko/stan/test_model.stan"
test_model = cmdstanpy.CmdStanModel(stan_file=test_model_path)

input_data = {
    "D": D,
    "N": N,
    "draws": len(posterior.coords["draw"]),
    "log_depths": np.log(test_tbl.sum(axis="sample")),
    "y": test_tbl.matrix_data.todense().T.astype(int),
    "x": dmat.values,
    "post_beta_var": post_beta_var,
    "post_phi": post_phi,
}

test_model_res = test_model.sample(
    data=input_data,
    iter_sampling=1,
    fixed_param=True
)
inf_test = az.from_cmdstanpy(
    test_model_res,
    dims={"all_log_lhood": ["post_draw", "tbl_sample", "class"]},
    coords={
        "post_draw": np.arange(len(post_beta_var)),
        "tbl_sample": test_tbl.ids(axis="sample"),
        "class": ["USA", "Malawi"]
    }
)

inf_test.to_netcdf("results/yatsuneko/inf.test.nc")
