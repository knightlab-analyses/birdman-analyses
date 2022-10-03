#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga/%x.%a.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --array=1-30

pwd; hostname; date

set -e

source ~/miniconda3/bin/activate birdman-analyses-final

echo Chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX

OUTDIR="/panfs/grahman/birdman-analyses-final/tcga/species/bacteria2"
LOGDIR="results/tcga/log/species/bacteria2"
TBL="data/tcga/processed/species/processed_tbl.bacteria.2.biom"
MD="data/tcga/processed/processed_md.2.tsv"
FORMULA="C(investigation, Treatment('TCGA-BRCA')) + race + gender"

mkdir -p $OUTDIR
mkdir -p $LOGDIR

echo Starting Python script...
time python src/tcga/run_birdman_chunked.py \
    --inference-dir $OUTDIR \
    --num-chunks $SLURM_ARRAY_TASK_MAX \
    --chunk-num $SLURM_ARRAY_TASK_ID \
    --table $TBL \
    --metadata $MD \
    --formula "${FORMULA}" \
    --chains 4 \
    --num-iter 500 \
    --num-warmup 1000 \
    --beta-prior 4.0 \
    --disp-scale 0.5 \
    --re-prior 3.0 \
    --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!
