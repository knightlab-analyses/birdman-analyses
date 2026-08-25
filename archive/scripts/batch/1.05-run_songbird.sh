#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/batch/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

source ~/miniconda3/bin/activate qiime2-2020.6

TABLE="results/batch/sim/sim_counts.biom"
MD="results/batch/sim/metadata.tsv"
RES_DIR="results/batch/songbird"
FORMULA="C(case_ctrl, Treatment(0))"

songbird multinomial \
    --input-biom $TABLE \
    --metadata-file $MD \
    --formula "${FORMULA}" \
    --random-seed 63 \
    --min-feature-count 0 \
    --min-sample-count 0 \
    --silent \
    --summary-dir $RES_DIR

