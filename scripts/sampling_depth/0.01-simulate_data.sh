#!/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/sampling_depth/%x.out
#SBATCH --partition=short
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00

pwd; hostname; date

set -e

source ~/miniconda3/bin/activate birdman-analyses-final

OUTDIR="results/sampling_depth/sim"
LOGDIR="results/sampling_depth/log"
mkdir -p $OUTDIR
mkdir -p $LOGDIR

NUM_SAMPS=500
MEAN_DEPTHS=(1000 2000 5000 10000 20000 50000 100000)
BETA_0=-8
BETA_1=3
LOGNORMAL_SD=0.2
INV_DISP=10

for D in ${MEAN_DEPTHS[@]}
do
    SIM_OUTPUT="${OUTDIR}/sim_counts.mean_depth_${D}.tsv"
    REG_OUTPUT="${OUTDIR}/reg_counts.mean_depth_${D}.nc"
    LOGFILE="${LOGDIR}/mean_depth_${D}.log"

    python src/sampling_depth/simulate_data.py \
        --sim-output $SIM_OUTPUT \
        --reg-output $REG_OUTPUT \
        --logfile $LOGFILE \
        --num-samples $NUM_SAMPS \
        --mean-depth $D \
        --lognormal-sd $LOGNORMAL_SD \
        --beta-0 $BETA_0 \
        --beta-1 $BETA_1 \
        --inv-disp $INV_DISP
done
