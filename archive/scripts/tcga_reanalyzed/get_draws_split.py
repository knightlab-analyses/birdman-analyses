#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=short
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G

import glob
from pathlib import Path
import re

import arviz as az
from joblib import delayed, Parallel

INF_PATH = Path("/panfs/grahman/birdman-analyses-final/tcga_reanalyzed/")
RES_PATH = Path("results/tcga_reanalyzed/split_draws")

invest_regex = re.compile("\[T\.TCGA-(\w+)\]")
feat_regex = re.compile("F\d{4}_(.*)\.nc")


def process_inf(inf_file, outdir):
    feat_name = feat_regex.search(inf_file).groups()[0]
    print(feat_name)
    post = az.from_netcdf(inf_file).posterior["beta_var"]
    post = (
        post
        .stack(sample=["chain", "draw"])
        .to_pandas()
        .T
        .filter(like="investigation")
    )
    post.columns = [invest_regex.search(x).groups()[0] for x in post.columns]
    post = post.reset_index(drop=True)

    outfile = outdir / f"{feat_name}.tsv"
    post.to_csv(outfile, sep="\t", index=True)


def get_draws(inf_dir, outdir):
    inf_file_list = sorted(list(map(str, inf_dir.glob("*.nc"))))

    (
        Parallel(n_jobs=20)
        (delayed(process_inf)(inf_file, outdir) for inf_file in inf_file_list)
    )


if __name__ == "__main__":
    for f in INF_PATH.iterdir():
        if "split" not in f.name:
            continue

        split_num = f.name[-1]
        res_dir = RES_PATH / f"split_{split_num}"
        res_dir.mkdir(parents=True, exist_ok=True)

        get_draws(f, res_dir)
