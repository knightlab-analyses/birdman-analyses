import os

import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
# repo root, so this runs from a fresh clone; override w/ BIRDMAN_ANALYSES_DIR
PROJ_DIR = os.environ.get(
    "BIRDMAN_ANALYSES_DIR",
    os.path.dirname(os.path.dirname(HERE))
)
MODEL_PATH = os.path.join(HERE, "stan", "abx_model_single.stan")
MD = pd.read_table(f"{PROJ_DIR}/data/relman_abx/processed/processed_md.tsv",
                   sep="\t", index_col=0)


class ABXModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        beta_prior: float = 10.0,
        disp_scale: float = 0.5,
        subj_prior: float = 3.0,
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

        # Levels of ABX treatment
        levels = [
            "preCp", "FirstCp", "FirstWPC", "Interim", "SecondCp",
            "SecondWPC", "PostCp"
        ]
        formula = f"C(antibiotic, Diff, levels={str(levels)})"
        self.create_regression(formula, MD)

        D = table.shape[0]
        A = np.log(1 / D)

        param_dict = {
            "depth": np.log(table.sum(axis="sample")),
            "num_subjs": len(self.subjects),
            "subj_map": samp_subj_map.values,
            "B_p": beta_prior,
            "disp_scale": disp_scale,
            "re_p": subj_prior,
            "A": A
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "phi", "subj_re"],
            dims={
                "beta_var": ["covariate"],
                "subj_re": ["subject", "covariate"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "covariate": self.colnames,
                "tbl_sample": self.sample_names,
                "subject": self.subjects
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"
        )
