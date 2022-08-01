from pkg_resources import resource_filename

import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

PROJ_DIR = "/home/grahman/projects/birdman-analyses-final"
MODEL_PATH = resource_filename("src", "obesity/stan/model_single.stan")
MD = pd.read_table(f"{PROJ_DIR}/data/obesity/processed/processed_md.tsv",
                   sep="\t", index_col=0)


class ModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        beta_prior: float = 3.0,
        disp_scale: float = 0.5,
        study_prior: float = 2.0,
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

        study_series = MD["study_id"].loc[self.sample_names]
        samp_study_map = study_series.astype("category").cat.codes + 1
        self.studies = np.sort(study_series.unique())

        formula = (
            "C(diet, Treatment('Standard')) + instrument"
        )
        self.create_regression(formula, MD)

        D = table.shape[0]
        A = np.log(1 / D)

        param_dict = {
            "depth": np.log(table.sum(axis="sample")),
            "num_studies": len(self.studies),
            "study_map": samp_study_map.values,
            "B_p": beta_prior,
            "disp_scale": disp_scale,
            "re_p": study_prior,
            "A": A
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "base_phi", "study_re", "study_disp"],
            dims={
                "beta_var": ["covariate"],
                "study_re": ["study"],
                "study_disp": ["study"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "covariate": self.colnames,
                "tbl_sample": self.sample_names,
                "study": self.studies
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"
        )
