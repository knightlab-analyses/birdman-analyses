#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/relman_abx/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

import glob
import re

import arviz as az
import pandas as pd

FEAT_REGEX = re.compile("F\d{4}_(.*).nc")

inference_dir = "/panfs/grahman/birdman-analyses-final/relman_abx/inferences"
outfile = "results/relman_abx/birdman_results.beta_var.tsv"
all_inf_files = glob.glob(f"{inference_dir}/*.nc")

def process_dataframe(df, feat_id, suffix=""):
    df = df.copy()
    df.reset_index(inplace=True, drop=True)
    df.columns.name = ""
    df.index = [feat_id]
    df.columns = [x + suffix for x in df.columns]
    return df

feat_diff_df_list = []
for inf_file in all_inf_files:
    this_feat_id = FEAT_REGEX.search(inf_file).groups()[0]
    this_feat_diff = az.from_netcdf(inf_file).posterior["beta_var"]
    this_feat_diff_mean = this_feat_diff.mean(["chain", "draw"]).to_dataframe().T
    this_feat_diff_std = this_feat_diff.std(["chain", "draw"]).to_dataframe().T

    this_feat_diff_mean = process_dataframe(
        this_feat_diff_mean,
        this_feat_id,
        suffix="_mean"
    )
    this_feat_diff_std = process_dataframe(
        this_feat_diff_std,
        this_feat_id,
        suffix="_std"
    )
    this_feat_diff_df = pd.concat(
        [this_feat_diff_mean, this_feat_diff_std],
        axis=1
    )
    feat_diff_df_list.append(this_feat_diff_df)

all_feat_diffs_df = pd.concat(feat_diff_df_list, axis=0)
all_feat_diffs_df.index.name = "Feature"
all_feat_diffs_df.to_csv(outfile, sep="\t", index=True)
