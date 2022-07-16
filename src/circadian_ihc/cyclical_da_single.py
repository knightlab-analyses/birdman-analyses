from pkg_resources import resource_filename

import biom
import cmdstanpy
import numpy as np
import pandas as pd

from birdman import SingleFeatureModel

PROJ_DIR = "/home/grahman/projects/birdman-analyses-final/"
MD = pd.read_table(f"{PROJ_DIR}/data/circadian_ihc/processed/circadianIHC_metadata_cleaned.txt",
                   sep="\t", index_col=0)
MODEL_PATH = resource_filename("src",
                               "circadian_ihc/stan/cyclical_da_single.stan")


class CircadianIHCModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        beta_prior: float = 2.0,
        inv_disp_sd: float = 0.5,
        subject_prior: float = 2.0,
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

        subj_series = MD["host_subject_id"].loc[self.sample_names]
        samp_subj_map = subj_series.astype("category").cat.codes + 1
        self.subjects = np.sort(subj_series.unique())

        self.create_regression("C(exposure_type, Treatment('Air'))", MD)

        D = table.shape[0]
        A = np.log(1 / D)

        param_dict = {
            "depth": np.log(table.sum(axis="sample")),
            "exposure": (MD["exposure_type"] == "IHC").astype(int).values,
            "S": len(self.subjects),
            "A": A,
            "time": MD["zt"].values,
            "subj_ids": samp_subj_map.values,
            "B_p": beta_prior,
            "inv_disp_sd": inv_disp_sd,
            "subj_prior": subject_prior,
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["reciprocal_phi", "beta_var", "subj_re",
                    "cos_coefs", "phase_shift", "period_offset_base",
                    "period_offset_ihc"],
            dims={
                "subj_re": ["subj_data", "subject"],
                "beta_var": ["covariate"],
                "cos_coefs": ["covariate"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "subject": self.subjects,
                "subj_data": ["intercept", "slope"],
                "covariate": self.colnames,
                "tbl_sample": self.sample_names
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"
        )
