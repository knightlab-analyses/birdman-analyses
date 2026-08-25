from pkg_resources import resource_filename

import arviz as az
import click
import cmdstanpy
import numpy as np
import pandas as pd

from src.logger import setup_loggers

rng = np.random.default_rng(63)


@click.command()
@click.option("--sim-output", required=True)
@click.option("--reg-output", required=True)
@click.option("--logfile", required=True)
@click.option("--num-samples", type=int)
@click.option("--mean-depth", type=int)
@click.option("--lognormal-sd", type=float)
@click.option("--beta-0", type=float)
@click.option("--beta-1", type=float)
@click.option("--inv-disp", type=float)
def simulate(
    sim_output,
    reg_output,
    logfile,
    num_samples,
    mean_depth,
    lognormal_sd,
    beta_0,
    beta_1,
    inv_disp,
):
    birdman_logger = setup_loggers(logfile)

    sim_model_path = resource_filename(
        "src", "sample_size/stan/sample_size_model_simulation_single.stan"
    )
    sim_model = cmdstanpy.CmdStanModel(stan_file=sim_model_path)

    mu = rng.lognormal(np.log(mean_depth), lognormal_sd, size=num_samples)
    log_depths = np.log(rng.poisson(mu))
    birdman_logger.info(f"Mean log-depth: {np.mean(log_depths)}")

    data = {
        "N": num_samples,
        "log_depths": log_depths,
        "beta_0": beta_0,
        "beta_1": beta_1,
        "inv_disp": inv_disp
    }
    birdman_logger.info(f"Sim data: {data}")

    sim_fit = sim_model.sample(
        data=data,
        iter_sampling=1,
        fixed_param=True,
        seed=63
    )

    sim_counts = sim_fit.stan_variable("sim_counts").squeeze()
    sim_counts = sim_counts.astype(int)

    columns = ["class_A", "class_B"]
    controls = [f"SA{x+1}" for x in range(int(num_samples/2))]
    cases = [f"SB{x+1}" for x in range(int(num_samples/2))]
    case_ctrl = [0] * int(num_samples/2) + [1] * int(num_samples/2)
    samples = controls + cases
    metadata = pd.DataFrame.from_dict({
        "count": sim_counts,
        "intercept": 1,
        "case_ctrl": case_ctrl,
        "log_depth": log_depths
    })
    metadata.index = samples
    birdman_logger.info(f"\n{metadata.head()}")
    metadata.to_csv(sim_output, sep="\t", index=True)
    x = metadata[["intercept", "case_ctrl"]].values

    reg_model_path = resource_filename(
        "src", "sample_size/stan/sample_size_model_regression_single.stan"
    )
    reg_model = cmdstanpy.CmdStanModel(stan_file=reg_model_path)

    reg_data = {
        "N": num_samples,
        "log_depths": log_depths,
        "counts": sim_counts,
        "x": x
    }

    reg = reg_model.sample(
        data=reg_data,
        iter_sampling=500,
        iter_warmup=1000,
        chains=4,
        seed=63
    )
    reg_inf = az.from_cmdstanpy(
        posterior=reg,
        dims={
            "beta_var": ["column"],
            "lam": ["tbl_sample"]
        },
        coords={
            "column": ["Intercept", "beta_case"],
            "tbl_sample": samples
        }
    )
    birdman_logger.info(f"{az.rhat(reg_inf)}")
    reg_inf.to_netcdf(reg_output)


if __name__ == "__main__":
    simulate()
