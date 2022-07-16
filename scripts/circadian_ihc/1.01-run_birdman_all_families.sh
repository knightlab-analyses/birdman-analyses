#!/usr/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/circadian_ihc/%x.out
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G

date; pwd; hostname

source ~/miniconda3/bin/activate birdman-analyses-final

time python src/circadian_ihc/run_cyclical_da_all_families.py
