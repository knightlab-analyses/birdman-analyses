#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/rare/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

pwd; hostname; date

set -e

source ~/miniconda3/bin/activate birdman-analyses-final

OUTDIR="results/rare/sim"
LOGDIR="results/rare/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR

NUM_SAMPS=300
NUM_FEATS=1000
BETA_0=-8
BETA_1=2

DEPTHS=(50 500 5000)
for d in "${DEPTHS[@]}"
do
    LOGFILE="${LOGDIR}/sim_${d}.log"

    python src/rare/simulate_data.py \
        --sim-output "${OUTDIR}/results_${d}.tsv" \
        --logfile $LOGFILE \
        --num-samples $NUM_SAMPS \
        --num-features $NUM_FEATS \
        --depth $d \
        --beta-0 $BETA_0 1 \
        --beta-1 $BETA_1 1
done

