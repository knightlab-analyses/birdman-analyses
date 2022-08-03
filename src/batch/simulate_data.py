from pkg_resources import resource_filename

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
@click.option("--reg-output", required=True)
@click.option("--logfile", required=True)
@click.option("--num-samples", type=int)
@click.option("--num-features", type=int)
@click.option("--mean-depth", type=int)
@click.option("--lognormal-depth-sd", type=float)
@click.option("--beta-0", type=(float, float))
@click.option("--beta-1", type=(float, float))
@click.option("--num-batches", type=int)
@click.option("--batch-offset-sd", type=float)
@click.option("--inv-disp", type=(float, float))
@click.option("--batch-disp-sd", type=float)
def simulate(
    sim_output,
    reg_output,
    logfile,
    num_samples,
    num_features,
    mean_depth,
    lognormal_depth_sd,
    beta_0,
    beta_1,
    num_batches,
    batch_offset_sd,
    inv_disp,
    batch_disp_sd
):
    birdman_logger = setup_loggers(logfile)

    sim_model_path = resource_filename(
        "src", "batch/stan/batch_model_simulation.stan"
    )
    sim_model = cmdstanpy.CmdStanModel(stan_file=sim_model_path)

    batches = np.arange(num_batches) + 1
    batch_offsets = rng.normal(0, batch_offset_sd,
                               size=(num_batches, num_features))
    batch_disps = rng.lognormal(0, batch_disp_sd,
                                size=(num_batches, num_features))
    batch_map = rng.choice(batches, size=num_samples)

    mu = rng.lognormal(np.log(mean_depth), lognormal_depth_sd, size=num_samples)
    log_depths = np.log(rng.poisson(mu))
    birdman_logger.info(f"Mean log-depth: {np.mean(log_depths)}")

    beta_0 = rng.normal(beta_0[0], beta_0[1], size=num_features)
    beta_1 = rng.normal(beta_1[0], beta_1[1], size=num_features)

    inv_disps = rng.lognormal(inv_disp[0], inv_disp[1], size=num_features)
    print(inv_disps)
    print(batch_disps)

    data = {
        "N": num_samples,
        "log_depths": log_depths,
        "D": num_features,
        "beta_0": beta_0,
        "beta_1": beta_1,
        "inv_disp": inv_disps,
        "num_batches": num_batches,
        "batch_map": batch_map,
        "batch_offsets": batch_offsets,
        "batch_disps": batch_disps
    }

    sim_fit = sim_model.sample(
        data=data,
        iter_sampling=1,
        fixed_param=True,
        seed=63
    )

    sim_counts = sim_fit.stan_variable("sim_counts").squeeze()
    sim_counts = sim_counts.astype(int)

    big_lam = sim_fit.stan_variable("big_lam").squeeze()
    big_alpha = sim_fit.stan_variable("big_alpha").squeeze()

    samp_ids = [f"S{x+1}" for x in range(num_samples)]
    feat_ids = [f"F{x+1}" for x in range(num_features)]
    tbl = biom.Table(
        sim_counts,
        sample_ids=samp_ids,
        observation_ids=feat_ids
    )

    big_alpha = pd.DataFrame(big_alpha)
    big_alpha.columns = samp_ids
    big_alpha.index = feat_ids
    big_alpha.to_csv("results/batch/sim/big_alpha.tsv", sep="\t", index=True)

    big_lam = pd.DataFrame(big_lam)
    big_lam.columns = samp_ids
    big_lam.index = feat_ids
    big_lam.to_csv("results/batch/sim/big_lam.tsv", sep="\t", index=True)

    # print(tbl.sum("sample"))
    # print(tbl.sum("observation"))

    with biom.util.biom_open(sim_output, "w") as f:
        tbl.to_hdf5(f, "sim")

    columns = ["class_A", "class_B"]
    controls = [f"SA{x+1}" for x in range(int(num_samples/2))]
    cases = [f"SB{x+1}" for x in range(int(num_samples/2))]
    case_ctrl = [0] * int(num_samples/2) + [1] * int(num_samples/2)
    samples = controls + cases
    metadata = pd.DataFrame.from_dict({
        "intercept": 1,
        "case_ctrl": case_ctrl,
        "log_depth": log_depths,
        "batch": batch_map
    })
    metadata.index = samples
    birdman_logger.info(f"\n{metadata.head()}")
    metadata.to_csv("results/batch/sim/metadata.tsv", sep="\t", index=True)

    x = metadata[["intercept", "case_ctrl"]].values

    # reg_model_path = resource_filename(
    #     "src", "sampling_depth/stan/depth_model_regression_single.stan"
    # )
    # reg_model = cmdstanpy.CmdStanModel(stan_file=reg_model_path)

    # reg_data = {
    #     "N": num_samples,
    #     "log_depths": log_depths,
    #     "counts": sim_counts,
    #     "x": x
    # }

    # reg = reg_model.sample(
    #     data=reg_data,
    #     iter_sampling=500,
    #     iter_warmup=1000,
    #     chains=4,
    #     seed=63
    # )
    # reg_inf = az.from_cmdstanpy(
    #     posterior=reg,
    #     dims={
    #         "beta_var": ["column"],
    #         "lam": ["tbl_sample"]
    #     },
    #     coords={
    #         "column": ["Intercept", "beta_case"],
    #         "tbl_sample": samples
    #     }
    # )
    # birdman_logger.info(f"{az.rhat(reg_inf)}")
    # reg_inf.to_netcdf(reg_output)


if __name__ == "__main__":
    simulate()
