#!/home/grahman/miniconda3/envs/birdman-benchmarking/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/relman_abx/%x.out
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

import biom
import numpy as np
import pandas as pd

tbl_file = "data/relman_abx/raw/45020_otu_table.biom"
tbl = biom.load_table(tbl_file)

md_file = "data/relman_abx/raw/494_20180101-113845.txt"
md = pd.read_table(md_file, sep="\t", index_col=0)

samps_in_common = set(md.index).intersection(tbl.ids())
tbl.filter(samps_in_common)
md = md.loc[samps_in_common]

# Filter observations by prevalence
n = 10
prevalence = tbl.pa(inplace=False).sum(axis="observation")
feats_to_keep = tbl.ids(axis="observation")[np.where(prevalence >= n)]
tbl.filter(ids_to_keep=feats_to_keep, axis="observation")

# Filter samples by depth
min_depth = 1000
depths = tbl.sum("sample")
samps_to_keep = tbl.ids()[np.where(depths >= min_depth)]
tbl.filter(samps_to_keep)
md = md.loc[samps_to_keep]

samp_map = {x: f"S{x}" for x in tbl.ids()}
tbl.update_ids(samp_map)
md.index = "S" + md.index

feat_map = {x: f"F{x}" for x in tbl.ids("observation")}
tbl.update_ids(feat_map, "observation")

# Save results
tbl_out_file = "data/relman_abx/processed/processed_tbl.biom"
with biom.util.biom_open(tbl_out_file, "w") as f:
    tbl.to_hdf5(f, "merge")

md.to_csv("data/relman_abx/processed/processed_md.tsv",
          sep="\t",
          index=True)
