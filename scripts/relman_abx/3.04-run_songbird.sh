#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/relman_abx/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00

source ~/miniconda3/bin/activate qiime2-2020.6

TABLE="data/relman_abx/processed/processed_tbl.biom"
MD="data/relman_abx/processed/processed_md.tsv"
RES_DIR="results/relman_abx/songbird"
FORMULA="host_subject_id + C(antibiotic, Diff, levels=['preCp', 'FirstCp', 'FirstWPC', 'Interim', 'SecondCp', 'SecondWPC', 'PostCp'])"

songbird multinomial \
    --input-biom $TABLE \
    --metadata-file $MD \
    --formula "${FORMULA}" \
    --random-seed 63 \
    --min-feature-count 0 \
    --min-sample-count 0 \
    --silent \
    --summary-dir $RES_DIR
