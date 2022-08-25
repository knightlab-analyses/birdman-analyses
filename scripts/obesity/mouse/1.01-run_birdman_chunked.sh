#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/mouse/%x.%a.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --array=1-20

pwd; hostname; date

set -e

source ~/miniconda3/bin/activate birdman-analyses-final

echo Chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX

OUTDIR="/panfs/grahman/birdman-analyses-final/obesity/mouse/inferences_genus"
LOGDIR="results/obesity/log/mouse"
TBL="data/obesity/processed/tbl_merged.genus.biom"
MD="data/obesity/processed/metadata.merged.tsv"
FORMULA="C(diet, Treatment('Standard')) + instrument"

mkdir -p $OUTDIR
mkdir -p $LOGDIR

echo Starting Python script...
time python src/obesity/run_birdman_chunked.py \
    --inference-dir $OUTDIR \
    --num-chunks $SLURM_ARRAY_TASK_MAX \
    --chunk-num $SLURM_ARRAY_TASK_ID \
    --table $TBL \
    --metadata $MD \
    --formula "${FORMULA}" \
    --chains 4 \
    --num-iter 500 \
    --num-warmup 1000 \
    --beta-prior 3.0 \
    --disp-scale 0.5 \
    --re-prior 3.0 \
    --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!
