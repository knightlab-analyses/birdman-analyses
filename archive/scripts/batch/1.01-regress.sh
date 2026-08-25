#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/batch/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

pwd; hostname; date

set -e

source ~/miniconda3/bin/activate birdman-analyses-final

OUTDIR="/panfs/grahman/birdman-analyses-final/batch/inferences"
LOGDIR="results/batch/log"

mkdir -p $OUTDIR
mkdir -p $LOGDIR

echo Starting Python script...
time python src/batch/regress_data.py \
    --inference-dir $OUTDIR \
    --chains 4 \
    --num-iter 500 \
    --num-warmup 1000 \
    --logfile "${LOGDIR}/regression.log" \
    --output "results/batch/regression_inf.nc" && echo Finished Python script!
