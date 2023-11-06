#!/home/grahman/miniconda3/envs/birdman-benchmarking/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga-re/%x.out
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

from functools import reduce
import biom
import numpy as np
import pandas as pd
from pathlib import Path
import re

tbl_dir = Path("data/tcga_reanalyzed/raw/")
blca_file = tbl_dir / "TableS8_BLCA.all.csv"
hnsc_file = tbl_dir / "TableS9_HNSC_all.csv"
brca_file = tbl_dir / "TableS10_BRCA_WGS.csv"

md_file = tbl_dir / "metadata_all_rep200_raw_15512samples.tsv"
md = pd.read_table(md_file, sep="\t", index_col=0)
md = md.query("experimental_strategy == 'WGS'")
md = md.set_index("knightlabID")
md = md[md["data_submitting_center_label"] != "MD Anderson - RPPA Core Facility (Proteomics)"]

md = md.query("gender != 'Not available'")
md = md.query("race != 'Not available'")

tbls = []
for f in [blca_file, hnsc_file, brca_file]:
    tbl = pd.read_table(f, sep=",", index_col=0)
    tbl = tbl.T
    tbl.index = [x.replace("g_", "") for x in tbl.index]
    samps_in_common = list(set(tbl.columns).intersection(md.index))
    tbl = tbl[samps_in_common]

    tbl = biom.Table(
        tbl.values,
        observation_ids=tbl.index,
        sample_ids=tbl.columns
    )
    tbls.append(tbl)

tbl_merged = reduce(lambda x, y: x.merge(y), tbls)

# Filter observations by prevalence
n = 20
prevalence = tbl_merged.pa(inplace=False).sum(axis="observation")
feats_to_keep = tbl_merged.ids(axis="observation")[np.where(prevalence >= n)]
tbl_merged.filter(ids_to_keep=feats_to_keep, axis="observation")

# Filter samples by depth
min_depth = 500
depths = tbl_merged.sum("sample")
samps_to_keep = tbl_merged.ids()[np.where(depths >= min_depth)]
tbl_merged.filter(samps_to_keep)
md = md.loc[samps_to_keep]

print(md.shape)
print(tbl_merged.shape)

# Save results
tbl_merged_out_file = "data/tcga_reanalyzed/processed/merged_tbl.biom"
with biom.util.biom_open(tbl_merged_out_file, "w") as f:
    tbl_merged.to_hdf5(f, "merge")

md.to_csv("data/tcga_reanalyzed/processed/processed_md.tsv",
          sep="\t",
          index=True)
