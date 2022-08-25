#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/mouse/%x.out
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

import os

import biom
import numpy as np
import pandas as pd


tbl_file = "data/obesity/raw/all.biom"
tbl = biom.load_table(tbl_file)
tbl_df = tbl.to_dataframe(dense=True)
tbl_df.index = "F" + tbl_df.index.astype(str)
tbl_df.columns = "S" + tbl_df.columns.astype(str)

md_file = "data/obesity/raw/all-metadata.txt"
md = pd.read_table(md_file, sep="\t", index_col=0)
md.index = "S" + md.index.astype(str)
md["study_id"] = "STUDY_" + md["study_id"].astype(str)

tax = pd.read_table("/databases/gg/13_8/taxonomy/97_otu_taxonomy.txt",
                    sep="\t", index_col=0, header=None)
tax.columns = ["Taxon"]
tax.index = "F" + tax.index.astype(str)
tax = tax["Taxon"].str.split("; ", expand=True)
tax.columns = list("kpcofgs")

tbl_df = tbl_df.join(tax, how="inner")
collapsed_df = tbl_df.groupby("g").sum()
collapsed_df = collapsed_df.drop(index=["g__"])

# Filter observations by prevalence
n = 10
prevalence = collapsed_df.clip(upper=1).sum(axis=1)
feats_to_keep = prevalence[prevalence >= n].index
collapsed_df = collapsed_df.loc[feats_to_keep]

# Filter samples by depth
d = 100
depths = collapsed_df.sum(axis=0)
samps_to_keep = depths[depths >= d].index
collapsed_df = collapsed_df.loc[:, samps_to_keep]

# Remove genera present in fewer than N studies
collapsed_df_joined = collapsed_df.T.join(md)
study_presence = (
    collapsed_df_joined
    .groupby("study_id")
    .sum()
    .clip(upper=1)
    .sum()
)
feats_to_keep = study_presence[study_presence >= 2].index
collapsed_df = collapsed_df.loc[feats_to_keep]

collapsed_tbl = biom.Table(
    collapsed_df.values,
    sample_ids=collapsed_df.columns,
    observation_ids=collapsed_df.index
)

md = md.loc[collapsed_df.columns]
na_samps = md[md["diet"].isna()]
diet_map = {True: "HFD", False: "Standard"}
md.loc[na_samps.index, "diet"] = [
    diet_map[x] for x in na_samps.index.str.contains("HFD")
]

# Save results
tbl_out_file = "data/obesity/processed/mouse/processed_tbl.genus.biom"
with biom.util.biom_open(tbl_out_file, "w") as f:
    collapsed_tbl.to_hdf5(f, "merge")

md.to_csv("data/obesity/processed/mouse/processed_md.tsv", sep="\t", index=True)

genus_counts = tbl_df["g"].value_counts().loc[collapsed_df.index]
genus_counts.name = "count"
genus_counts.to_csv("results/obesity/mouse/genus_counts.tsv", sep="\t", index=True)

study_dir = "data/obesity/processed/study_tables"
os.makedirs(study_dir, exist_ok=True)
for study, study_md in md.groupby("study_id"):
    idx = study_md.index
    tbl_study = collapsed_tbl.filter(idx, inplace=False)
    tbl_study.remove_empty()
    outpath = os.path.join(study_dir, f"{study}_tbl.genus.biom")

    with biom.util.biom_open(outpath, "w") as f:
        tbl_study.to_hdf5(f, "filt")

    print(f"{study}: {tbl_study.shape}")
