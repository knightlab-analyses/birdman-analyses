import os

import arviz as az
import biom
import click
import cmdstanpy
import numpy as np
import pandas as pd

from src.logger import setup_loggers

rng = np.random.default_rng(63)

@click.command()
@click.option("--sim-output", required=True)
@click.option("--logfile", required=True)
@click.option("--num-samples", type=int)
@click.option("--num-features", type=int)
@click.option("--depth", type=int)
@click.option("--beta-0", type=(float, float))
@click.option("--beta-1", type=(float, float))
def simulate(
    sim_output,
    logfile,
    num_samples,
    num_features,
    depth,
    beta_0,
    beta_1,
):
    birdman_logger = setup_loggers(logfile)

    beta_0_truth = rng.normal(beta_0[0], beta_0[1], size=num_features)
    beta_1_truth = rng.normal(beta_1[0], beta_1[1], size=num_features)
    birdman_logger.info(f"beta_0: {beta_0_truth}")
    birdman_logger.info(f"beta_1: {beta_1_truth}")

    model_path = "src/rare/stan/simulate.stan"
    model = cmdstanpy.CmdStanModel(stan_file=model_path)
    data = {
        "N": num_samples,
        "D": num_features,
        "depth": depth,
        "beta_0": beta_0_truth,
        "beta_1": beta_1_truth
    }

    sim_fit = model.sample(
        data=data,
        iter_sampling=1,
        fixed_param=True,
        seed=63,
    )
    sim_counts = sim_fit.stan_variable("sim_counts").squeeze()

    reg_model_path = "src/rare/stan/regress_single.stan"
    reg_model = cmdstanpy.CmdStanModel(stan_file=reg_model_path)

    md = np.concatenate([np.zeros(num_samples//2), np.ones(num_samples//2)])
    all_vals = []

    for i in range(num_features):
        data = {
            "N": num_samples,
            "y": sim_counts[:, i].astype(int),
            "x": md,
            "beta_0": beta_0_truth[i],
            "depth": depth
        }
        reg_fit = reg_model.sample(
            data=data,
            iter_warmup=500,
            iter_sampling=500,
            seed=63
        )
        beta_1_fit = reg_fit.stan_variable("beta_1")
        beta_1_mean, beta_1_var = beta_1_fit.mean(), beta_1_fit.var()
        vals = pd.Series([beta_1_truth[i], beta_1_mean, beta_1_var, depth])
        all_vals.append(vals)

    val_df = pd.DataFrame(all_vals)
    val_df.columns = [
        "truth",
        "est_mean",
        "est_var",
        "depth"
    ]
    val_df.to_csv(sim_output, sep="\t")


if __name__ == "__main__":
    simulate()
