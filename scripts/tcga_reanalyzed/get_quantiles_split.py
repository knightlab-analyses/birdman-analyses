#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=short
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G

from pathlib import Path

import pandas as pd
from joblib import delayed, Parallel


def process_feat(feat_draws_file):
    feat_name = feat_draws_file.stem
    feat_df = pd.read_table(feat_draws_file, sep="\t", index_col=0)

    quantiles = feat_df.quantile([0.025, 0.975]).reset_index()
    quantiles["feature"] = feat_name
    quantiles = quantiles.rename(columns={"index": "quantile"})
    return quantiles


if __name__ == "__main__":
    draw_path = Path("results/tcga_reanalyzed/split_draws")
    for split_dir in draw_path.iterdir():
        split_num = split_dir.name[-1]
        draws_files = split_dir.glob("*.tsv")

        feat_dfs = Parallel(n_jobs=20)(delayed(process_feat)(f) for f in draws_files)
        quantile_df = pd.concat(feat_dfs)

        outfile = f"results/tcga_reanalyzed/feature_quantiles.split_{split_num}.tsv"
        quantile_df.to_csv(outfile, sep="\t", index=False)
