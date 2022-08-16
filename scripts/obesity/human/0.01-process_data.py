#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/human/%x.out
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

import functools
import pathlib

import biom
import numpy as np
import pandas as pd

md_12879_path = "data/obesity/raw/qiita_12879/12879_20210429-100954.txt"
na_values = ["not applicable", "not provided"]
md_12879 = pd.read_table(md_12879_path, sep="\t", index_col=0,
                         na_values=na_values).dropna()
md_12879.index = [f"S{x}" for x in md_12879.index]

def get_first_timept(x):
    return x.sort_values(by="collection_timestamp").head(1)

md_12879_first = (
    md_12879
    .groupby("host_subject_id")
    .apply(get_first_timept)
    .droplevel(0)
    .assign(obese=lambda x:
        (x["host_body_mass_index"] > 35).map({True: "Obese", False: "Lean"})
    )
)

study_12879_path = pathlib.Path("data/obesity/raw/qiita_12879")
study_12879_tbls = [biom.load_table(f) for f in study_12879_path.glob("*.biom")]

tbl_12879 = functools.reduce(lambda a, b: a.merge(b), study_12879_tbls)
tbl_12879.update_ids({x: f"S{x}" for x in tbl_12879.ids()})
samps_to_keep = list(set(tbl_12879.ids()).intersection(md_12879_first.index))
md_12879_first = md_12879_first.loc[samps_to_keep]
tbl_12879.filter(samps_to_keep)
tbl_12879.remove_empty()

md_77_path = "data/obesity/raw/qiita_77/77_20180101-113930.txt"
md_77 = pd.read_table(md_77_path, sep="\t", index_col=0)
md_77 = (
    md_77
    .query("obesitycat == ['Obese', 'Lean']")
    .rename(columns={"obesitycat": "obese"})
)
md_77.index = [f"S{x}" for x in md_77.index]

def get_first_sampledate(x):
    return x.query("sampledate == 'TimePoint1'")

md_77_first = (
    md_77
    .groupby("host_subject_id")
    .apply(get_first_sampledate)
    .droplevel(0)
)

tbl_77 = biom.load_table("data/obesity/raw/qiita_77/45469_otu_table.biom")
tbl_77.update_ids({x: f"S{x}" for x in tbl_77.ids()})
samps_to_keep = list(set(md_77_first.index).intersection(tbl_77.ids()))
tbl_77.filter(samps_to_keep)
md_77 = md_77.loc[samps_to_keep]
tbl_77.remove_empty()

cols = ["obese", "qiita_study_id"]
concat_md = pd.concat([md_77[cols], md_12879_first[cols]])
concat_tbl = tbl_12879.merge(tbl_77, observation="intersection")
concat_tbl.remove_empty()
concat_tbl_df = concat_tbl.to_dataframe(dense=True)
concat_tbl_df.index = "F" + concat_tbl_df.index.astype(str)

tax = pd.read_table("/databases/gg/13_8/taxonomy/97_otu_taxonomy.txt",
                    sep="\t", index_col=0, header=None)
tax.columns = ["Taxon"]
tax.index = "F" + tax.index.astype(str)
tax = tax["Taxon"].str.split("; ", expand=True)
tax.columns = list("kpcofgs")

concat_tbl_df = concat_tbl_df.join(tax, how="inner")
collapsed_df = concat_tbl_df.groupby("g").sum()
collapsed_df = collapsed_df.drop(index=["g__"])

mouse_tbl_path = "data/obesity/processed/processed_tbl.genus.biom"
genera_to_keep = list(set(
    biom.load_table(mouse_tbl_path)
    .ids("observation")
).intersection(collapsed_df.index))
collapsed_df = collapsed_df.loc[genera_to_keep]

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

collapsed_tbl = biom.Table(
    collapsed_df.values,
    sample_ids=collapsed_df.columns,
    observation_ids=collapsed_df.index
)

tbl_out_file = "data/obesity/processed/human/processed_tbl.genus.biom"
with biom.util.biom_open(tbl_out_file, "w") as f:
    collapsed_tbl.to_hdf5(f, "filtered")

concat_md = concat_md.loc[collapsed_tbl.ids()]
md_out_file = "data/obesity/processed/human/processed_md.tsv"
concat_md.to_csv(md_out_file, sep="\t", index=True)

genus_counts = concat_tbl_df["g"].value_counts().loc[collapsed_df.index]
genus_counts.name = "count"
genus_counts.to_csv("results/obesity/genus_counts.human.tsv", sep="\t", index=True)
