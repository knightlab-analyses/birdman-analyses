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

tbl1 = biom.load_table("data/obesity/processed/mouse/processed_tbl.genus.biom")
md1 = pd.read_table("data/obesity/processed/mouse/processed_md.tsv", sep="\t",
                    index_col=0)

tbl2 = biom.load_table("data/obesity/intermediate/study_11548.biom")
tbl2_df = tbl2.to_dataframe(dense=True)
tbl2_df.index = "F" + tbl2_df.index.astype(str)
tbl2_df.columns = "S" + tbl2_df.columns.astype(str)
md2 = pd.read_table("data/obesity/intermediate/study_11548.tsv", sep="\t",
                    index_col=0)

def get_first_timepoint(df):
    return df.sort_values(by="collection_timestamp").head(1)

cols_to_keep = ["host_subject_id", "diet", "variable_region", "instrument",
                "study_id"]
first_timepoint = (
    md2
    .query("exposure_type == 'Air'")
    .groupby("host_subject_id")
    .apply(get_first_timepoint)
    .reset_index(level=1)
    .set_index("#SampleID")
    .assign(
        diet=lambda x: x["diet"].replace({"Regular chow": "Standard"}),
        study_id="STUDY_11548",
        instrument="Illumina",
        variable_region="V4"
    )
    [cols_to_keep]
)
first_timepoint.index = [f"S{x}" for x in first_timepoint.index]
tbl2_df = tbl2_df[first_timepoint.index]

tax = pd.read_table("/databases/gg/13_8/taxonomy/97_otu_taxonomy.txt",
                    sep="\t", index_col=0, header=None)
tax.columns = ["Taxon"]
tax.index = "F" + tax.index.astype(str)
tax = tax["Taxon"].str.split("; ", expand=True)
tax.columns = list("kpcofgs")

tbl2_df = tbl2_df.join(tax, how="inner")
collapsed_df = tbl2_df.groupby("g").sum()
collapsed_df = collapsed_df.drop(index=["g__"])

samps_to_keep = md1.query("study_id != 'STUDY_11548'").index
tbl1.filter(samps_to_keep)

collapsed_tbl = biom.Table(
    collapsed_df.values,
    sample_ids=collapsed_df.columns,
    observation_ids=collapsed_df.index
)

tbl_merged = tbl1.merge(collapsed_tbl, observation="intersection")
md_merged = pd.concat([md1.loc[samps_to_keep], first_timepoint])

tbl_outfile = "data/obesity/processed/mouse/tbl_merged.genus.biom"
with biom.util.biom_open(tbl_outfile, "w") as f:
    tbl_merged.to_hdf5(f, "merged")

md_outfile = "data/obesity/processed/mouse/metadata.merged.tsv"
md_merged.to_csv(md_outfile, sep="\t", index=True)
