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

OUTDIR="results/batch/sim"
LOGDIR="results/batch/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR

NUM_SAMPS=300
NUM_FEATS=10
MEAN_DEPTH=5000
BETA_0=-8
BETA_1=2
NUM_BATCHES=10
LOGNORMAL_SD=0.2

LOGFILE="${LOGDIR}/sim.log"

python src/batch/simulate_data.py \
    --sim-output $OUTDIR \
    --logfile $LOGFILE \
    --num-samples $NUM_SAMPS \
    --num-features $NUM_FEATS \
    --mean-depth $MEAN_DEPTH \
    --lognormal-depth-sd $LOGNORMAL_SD \
    --beta-0 $BETA_0 1 \
    --beta-1 $BETA_1 1 \
    --num-batches $NUM_BATCHES
