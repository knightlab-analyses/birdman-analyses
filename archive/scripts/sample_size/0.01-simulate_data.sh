#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/sample_size/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

pwd; hostname; date

set -e

source ~/miniconda3/bin/activate birdman-analyses-final

OUTDIR="results/sample_size/sim"
LOGDIR="results/sample_size/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR

NUM_SAMPS=(50 100 200 500 1000 5000 10000)
MEAN_DEPTH=50000
BETA_0=-8
BETA_1=3
LOGNORMAL_SD=0.5
INV_DISP=10

for N in ${NUM_SAMPS[@]}
do
    SIM_OUTPUT="${OUTDIR}/sim_counts.num_samps_${N}.tsv"
    REG_OUTPUT="${OUTDIR}/reg_counts.num_samps_${N}.nc"
    LOGFILE="${LOGDIR}/num_samps_${N}.log"

    python src/sample_size/simulate_data.py \
        --sim-output $SIM_OUTPUT \
        --reg-output $REG_OUTPUT \
        --logfile $LOGFILE \
        --num-samples $N \
        --mean-depth $MEAN_DEPTH \
        --lognormal-sd $LOGNORMAL_SD \
        --beta-0 $BETA_0 \
        --beta-1 $BETA_1 \
        --inv-disp $INV_DISP
done
