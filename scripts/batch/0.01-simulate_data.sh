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

NUM_SAMPS=100
NUM_FEATS=10
MEAN_DEPTH=5000
BETA_0=-8
BETA_1=2
NUM_BATCHES=10
LOGNORMAL_SD=0.2
INV_DISP=-2
INV_DISP_SD=0.5
BATCH_DISP_SD=0.2

SIM_OUTPUT="${OUTDIR}/sim_counts.biom"
REG_OUTPUT="${OUTDIR}/reg_counts.nc"
LOGFILE="${LOGDIR}/sim.log"

python src/batch/simulate_data.py \
    --sim-output $SIM_OUTPUT \
    --reg-output $REG_OUTPUT \
    --logfile $LOGFILE \
    --num-samples $NUM_SAMPS \
    --num-features $NUM_FEATS \
    --mean-depth $MEAN_DEPTH \
    --lognormal-depth-sd $LOGNORMAL_SD \
    --beta-0 $BETA_0 1 \
    --beta-1 $BETA_1 1 \
    --num-batches $NUM_BATCHES \
    --batch-offset-sd 2 \
    --inv-disp $INV_DISP $INV_DISP_SD \
    --batch-disp-sd $BATCH_DISP_SD
