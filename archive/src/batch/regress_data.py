import os
import json
from tempfile import TemporaryDirectory
import time

import arviz as az
import biom
from birdman import ModelIterator
from birdman.model_util import concatenate_inferences
import cmdstanpy
import click
import numpy as np
import pandas as pd

from src.logger import setup_loggers
from src.batch.regression_model_single import RegressionModelSingle

PROJ_DIR = "/home/grahman/projects/birdman-analyses-final"
TABLE_FILE = f"{PROJ_DIR}/results/batch/sim/sim_counts.biom"
TABLE = biom.load_table(TABLE_FILE)
FIDS = TABLE.ids(axis="observation")


@click.command()
@click.option("--inference-dir", required=True)
@click.option("--chains", default=4)
@click.option("--num-iter", default=500)
@click.option("--num-warmup", default=1000)
@click.option("--logfile", required=True)
@click.option("--output", required=True)
def run_birdman(
    inference_dir,
    chains,
    num_iter,
    num_warmup,
    logfile,
    output
):
    birdman_logger = setup_loggers(logfile)

    model_iter = ModelIterator(
        TABLE,
        RegressionModelSingle,
        chains=chains,
        num_iter=num_iter,
        num_warmup=num_warmup,
    )

    infs = []
    for feature_id, model in model_iter:
        feature_num = np.where(FIDS == feature_id)[0].item()
        feature_num_str = str(feature_num).zfill(4)
        birdman_logger.info(f"Feature num: {feature_num_str}")
        birdman_logger.info(f"Feature ID: {feature_id}")

        tmpdir = f"{inference_dir}/tmp/F{feature_num_str}_{feature_id}"

        os.makedirs(tmpdir, exist_ok=True)

        with TemporaryDirectory(dir=tmpdir) as t:
            model.compile_model()
            model.fit_model(sampler_args={"output_dir": t})

            inf = model.to_inference()
            birdman_logger.info(inf.posterior)

            rhat = az.rhat(inf)
            birdman_logger.info("Rhat:")
            birdman_logger.info(rhat)
            if (rhat > 1.05).to_array().any().item():
                birdman_logger.warning(
                    f"{feature_id} has Rhat values > 1.05"
                )

            birdman_logger.info(f"Finished fitting!")
            time.sleep(10)

        infs.append(inf)

    concat_inf = concatenate_inferences(
        infs,
        coords={"feature": TABLE.ids("observation")}
    )
    concat_inf.to_netcdf(output)


if __name__ == "__main__":
    run_birdman()

