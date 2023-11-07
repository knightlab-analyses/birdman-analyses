#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga/%x.out
#SBATCH --partition=long
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00

source ~/miniconda3/bin/activate qiime2-2020.6

TABLE="data/tcga/processed/species/processed_tbl.bacteria.2.biom"
MD="data/tcga/processed/processed_md.2.tsv"
RES_DIR="results/tcga/songbird"
FORMULA="C(investigation, Treatment('TCGA-BRCA')) + race + gender"

songbird multinomial \
    --input-biom $TABLE \
    --metadata-file $MD \
    --formula "${FORMULA}" \
    --random-seed 63 \
    --min-feature-count 0 \
    --min-sample-count 0 \
    --silent \
    --summary-dir $RES_DIR
