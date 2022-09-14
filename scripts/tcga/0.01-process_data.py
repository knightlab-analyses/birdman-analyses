#!/home/grahman/miniconda3/envs/birdman-benchmarking/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga/%x.out
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

import biom
import numpy as np
import pandas as pd

tbl_file = "data/tcga/raw/species/immune_rep200_counts_bacteria_all_primary_tumors.biom"
tbl = biom.load_table(tbl_file)

md_file = "data/tcga/raw/metadata_immune_WGS_AllSeqCenters_Primary_Tumor.txt"
md = pd.read_table(md_file, sep="\t", index_col=0,
                   na_values=["Not available", "None"])
md = md.drop_duplicates(subset=["tcga_sample_id"], keep="first")

md = md[md["data_submitting_center_label"] != "MD Anderson - RPPA Core Facility (Proteomics)"]
md = md[md["race"].isin(["WHITE", "BLACK OR AFRICAN AMERICAN", "ASIAN"])]
md = md[md["vital_status"] == "Alive"]
md = md.dropna(subset=[
    "race",
    "investigation",
    "data_submitting_center_label",
    "tissue_source_site_label"
])

samps_in_common = set(md.index).intersection(tbl.ids())
tbl.filter(samps_in_common)
md = md.loc[samps_in_common]

# Filter observations by prevalence
n = 100
prevalence = tbl.pa(inplace=False).sum(axis="observation")
feats_to_keep = tbl.ids(axis="observation")[np.where(prevalence >= n)]
tbl.filter(ids_to_keep=feats_to_keep, axis="observation")

# Filter samples by depth
min_depth = 500
depths = tbl.sum("sample")
samps_to_keep = tbl.ids()[np.where(depths >= min_depth)]
tbl.filter(samps_to_keep)
md = md.loc[samps_to_keep]

samp_map = {x: f"S{x}" for x in tbl.ids()}
tbl.update_ids(samp_map)
md.index = "S" + md.index

invest_vc = md["investigation"].value_counts()
invest_to_keep = invest_vc[invest_vc >= 20].index

md = md[md["investigation"].isin(invest_to_keep)]
tbl.filter(md.index)

# Save results
tbl_out_file = "data/tcga/processed/species/processed_tbl.bacteria.biom"
with biom.util.biom_open(tbl_out_file, "w") as f:
    tbl.to_hdf5(f, "merge")

md.to_csv("data/tcga/processed/processed_md.tsv",
          sep="\t",
          index=True)
