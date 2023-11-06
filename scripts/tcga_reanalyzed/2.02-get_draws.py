#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=short
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G

import glob
import re

import arviz as az
from joblib import delayed, Parallel


inf_dir = "/panfs/grahman/birdman-analyses-final/tcga_reanalyzed/v1"
inf_file_list = sorted(glob.glob(f"{inf_dir}/*.nc"))

post_list = []

invest_regex = re.compile("\[T\.TCGA-(\w+)\]")
feat_regex = re.compile("F\d{4}_(.*)\.nc")


def process_inf(inf_file):
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

    post.to_csv(
        f"results/tcga_reanalyzed/draws/{feat_name}.tsv",
        sep="\t",
        index=True
    )

(
    Parallel(n_jobs=20)
    (delayed(process_inf)(inf_file) for inf_file in inf_file_list)
)
