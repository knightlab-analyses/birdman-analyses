#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.out
#SBATCH --partition=long
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00

source ~/miniconda3/bin/activate qiime2-2020.6

MD="data/tcga_reanalyzed/processed/processed_md.tsv"
RES_DIR="results/tcga_reanalyzed/songbird_split"
FORMULA="C(investigation, Treatment('TCGA-BRCA')) + race + gender + data_submitting_center_label"

for i in 0 1 2 3 4
do
    TBL="data/tcga_reanalyzed/processed/split_${i}/merged_tbl.train.${i}.biom"
    OUTDIR="${RES_DIR}/split_${i}"

    mkdir $OUTDIR

    songbird multinomial \
        --input-biom $TBL \
        --metadata-file $MD \
        --formula "${FORMULA}" \
        --random-seed 42 \
        --differential-prior 4.0 \
        --epochs 1000 \
        --min-feature-count 0 \
        --min-sample-count 0 \
        --silent \
        --summary-dir $OUTDIR
done
