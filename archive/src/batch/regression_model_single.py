import json
from pkg_resources import resource_filename

import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

PROJ_DIR = "/home/grahman/projects/birdman-analyses-final"
MODEL_PATH = resource_filename("src", "batch/stan/batch_model_regression.stan")
MD = pd.read_table(f"{PROJ_DIR}/results/batch/sim/metadata.tsv",
                   sep="\t", index_col=0)
PARAMS_FILE= f"{PROJ_DIR}/results/batch/sim/params.json"
with open(PARAMS_FILE, "r") as f:
    PARAMS = json.load(f)

class RegressionModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        num_iter: int,
        num_warmup: int,
        **kwargs
    ):
        super().__init__(
            table=table,
            feature_id=feature_id,
            model_path=MODEL_PATH,
            num_iter=num_iter,
            num_warmup=num_warmup,
            **kwargs
        )

        batch_series = MD["batch"].loc[self.sample_names]
        samp_batch_map = batch_series.astype("category").cat.codes + 1
        self.batches = np.sort(batch_series.unique())

        formula = f"C(case_ctrl, Treatment(0))"
        self.create_regression(formula, MD)

        param_dict = {
            "log_depths": PARAMS["log_depths"],
            "num_batches": len(self.batches),
            "batch_map": samp_batch_map.values,
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "inv_disp", "batch_offsets", "batch_disps"],
            dims={
                "beta_var": ["covariate"],
                "batch_offsets": ["batch"],
                "batch_disps": ["batch"],
            },
            coords={
                "covariate": self.colnames,
                "batch": self.batches
            },
        )
