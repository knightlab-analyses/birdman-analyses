from pkg_resources import resource_filename

import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

PROJ_DIR = "/home/grahman/projects/birdman-analyses-final"
MODEL_PATH = resource_filename("src", "tcga/stan/model_single.stan")


class ModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        metadata: pd.DataFrame,
        formula: str,
        beta_prior: float = 3.0,
        disp_scale: float = 0.5,
        re_prior: float = 2.0,
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

        MD = metadata

        center_series = MD["data_submitting_center_label"].loc[self.sample_names]
        samp_center_map = center_series.astype("category").cat.codes + 1
        self.centers = np.sort(center_series.unique())

        self.create_regression(formula, MD)

        D = table.shape[0]
        A = np.log(1 / D)

        param_dict = {
            "depth": np.log(table.sum(axis="sample")),
            "num_center": len(self.centers),
            "center_map": samp_center_map.values,
            "B_p": beta_prior,
            "disp_scale": disp_scale,
            "re_p": re_prior,
            "A": A
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "base_phi", "center_re", "center_disp"],
            dims={
                "beta_var": ["covariate"],
                "center_re": ["center"],
                "center_disp": ["center"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "covariate": self.colnames,
                "tbl_sample": self.sample_names,
                "center": self.centers
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"
        )
