from tempfile import TemporaryDirectory
import time

import arviz as az
import biom
from birdman.model_util import concatenate_inferences

from src.circadian_ihc.cyclical_da_single import CircadianIHCModelSingle

table = biom.load_table("data/circadian_ihc/processed/tbl.family.filt.biom")
infs = []

for family in table.ids(axis="observation"):
    with TemporaryDirectory(dir="/panfs/grahman/tmp") as tmpdir:
        model = CircadianIHCModelSingle(
            table,
            feature_id=family,
            seed=63,
            beta_prior = 2.0,
            inv_disp_sd = 0.5,
            subject_prior = 2.0,
            num_iter = 500,
            num_warmup = 1000,
        )
        model.compile_model()
        model.fit_model(sampler_args={"output_dir": tmpdir})
        inf = model.to_inference_object()
        infs.append(inf)

    time.sleep(5)

all_families_inf = concatenate_inferences(
    infs,
    {"family": table.ids(axis="observation")},
    "family"
)
all_families_inf.to_netcdf("results/circadian_ihc/all_families_inf.nc")
