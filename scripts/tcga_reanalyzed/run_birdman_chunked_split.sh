#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/tcga_reanalyzed/%x.%a.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --array=1-20

pwd; hostname; date

set -e

source ~/miniconda3/bin/activate birdman-analyses-final

echo Chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX

for i in 0 1 2 3 4
do
    OUTDIR="/panfs/grahman/birdman-analyses-final/tcga_reanalyzed/split_${i}"
    LOGDIR="results/tcga_reanalyzed/log/split_${i}"
    TBL="data/tcga_reanalyzed/processed/split_${i}/merged_tbl.train.${i}.biom"
    MD="data/tcga_reanalyzed/processed/processed_md.tsv"
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
        --num-iter 300 \
        --num-warmup 1000 \
        --beta-prior 4.0 \
        --disp-scale 0.5 \
        --re-prior 3.0 \
        --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!
done
