#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=long
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00

source ~/miniconda3/bin/activate qiime2-2020.6

TABLE="data/tcga_reanalyzed/processed/merged_tbl.biom"
MD="data/tcga_reanalyzed/processed/processed_md.tsv"
RES_DIR="results/tcga_reanalyzed/songbird"
FORMULA="C(investigation, Treatment('TCGA-BRCA')) + race + gender + data_submitting_center_label"

songbird multinomial \
    --input-biom $TABLE \
    --metadata-file $MD \
    --formula "${FORMULA}" \
    --random-seed 42 \
    --differential-prior 4.0 \
    --epochs 1000 \
    --min-feature-count 0 \
    --min-sample-count 0 \
    --silent \
    --summary-dir $RES_DIR
