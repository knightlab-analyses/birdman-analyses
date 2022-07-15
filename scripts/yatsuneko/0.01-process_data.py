#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/yatsuneko/%x.out
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

import biom
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


tbl_file = "data/yatsuneko/raw/134451_reference-hit.biom"
tbl = biom.load_table(tbl_file)

md_file = "data/yatsuneko/raw/850_20220314-110911.txt"
md = pd.read_table(md_file, sep="\t", index_col=0)
countries = ["Malawi", "USA"]
md = md.query("geo_loc_name == @countries")

samps_in_common = set(md.index).intersection(tbl.ids())
tbl.filter(samps_in_common)
md = md.loc[samps_in_common]

# Filter observations by prevalence
n = 50
prevalence = tbl.pa(inplace=False).sum(axis="observation")
feats_to_keep = tbl.ids(axis="observation")[np.where(prevalence >= n)]
tbl.filter(ids_to_keep=feats_to_keep, axis="observation")

samp_map = {x: f"S{x}" for x in tbl.ids()}
tbl.update_ids(samp_map)
md.index = "S" + md.index

# Save results
tbl_out_file = "data/yatsuneko/processed/processed_tbl.biom"
with biom.util.biom_open(tbl_out_file, "w") as f:
    tbl.to_hdf5(f, "merge")

md.to_csv("data/yatsuneko/processed/processed_md.tsv",
          sep="\t",
          index=True)

# Train-Test
# Need to make sure sample order is retained
train_md, test_md = train_test_split(md, test_size=0.3, random_state=42)

train_tbl = tbl.filter(train_md.index, inplace=False)
test_tbl = tbl.filter(test_md.index, inplace=False)

train_md = train_md.loc[train_tbl.ids()]
test_md = test_md.loc[test_tbl.ids()]

train_md.to_csv("data/yatsuneko/processed/processed_md.train.tsv",
                sep="\t",
                index=True)
test_md.to_csv("data/yatsuneko/processed/processed_md.test.tsv",
            sep="\t",
            index=True)

train_tbl_out_file = "data/yatsuneko/processed/processed_tbl.train.biom"
with biom.util.biom_open(train_tbl_out_file, "w") as f:
    train_tbl.to_hdf5(f, "merge")

test_tbl_out_file = "data/yatsuneko/processed/processed_tbl.test.biom"
with biom.util.biom_open(test_tbl_out_file, "w") as f:
    test_tbl.to_hdf5(f, "merge")
