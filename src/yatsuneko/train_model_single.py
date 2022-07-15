from pkg_resources import resource_filename

import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

PROJ_DIR = "/home/grahman/projects/birdman-analyses-final"
MODEL_PATH = resource_filename("src", "yatsuneko/stan/train_model.stan")
MD = pd.read_table(f"{PROJ_DIR}/data/yatsuneko/processed/processed_md.tsv",
                   sep="\t", index_col=0)


class TrainModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        disp_scale: float = 0.5,
        beta_prior: float = 2.0,
        num_iter: int = 500,
        num_warmup: int = 1000,
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

        self.create_regression("C(geo_loc_name, Treatment('USA'))", MD)
        D = table.shape[0]
        A = np.log(1 / D)

        param_dict = {
            "log_depth": np.log(table.sum(axis="sample")),
            "disp_scale": disp_scale,
            "B_p": beta_prior,
            "A": A,
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "phi"],
            dims={
                "beta_var": ["covariate"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"],
            },
            coords={
                "covariate": self.colnames,
                "tbl_sample": self.sample_names,
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"
        )
