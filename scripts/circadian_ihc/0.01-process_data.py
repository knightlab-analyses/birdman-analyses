#!/home/grahman/miniconda3/envs/birdman-analyses-final/bin/python
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/circadian_ihc/%x.out
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

import biom
import numpy as np
import pandas as pd

tbl_file = "data/circadian_ihc/raw/CircadianIHC_family_counts.biom"
md_file = "data/circadian_ihc/processed/circadianIHC_metadata_cleaned.txt"

md = pd.read_table(md_file, sep="\t", index_col=0)

tbl = (
    biom.load_table(tbl_file)
    .filter(md.index, axis="sample", inplace=False)
)

p = 30
feat_prev = tbl.pa(inplace=False).sum(axis="observation")
feats_to_keep = tbl.ids(axis="observation")[np.where(feat_prev >= p)]
tbl.filter([x for x in feats_to_keep if "f__" in x],
            axis="observation")

out_tbl_file = "data/processed/tbl.family.filt.biom"

with biom.util.biom_open(out_tbl_file, "w") as f:
    tbl.to_hdf5(f, "filtered")
