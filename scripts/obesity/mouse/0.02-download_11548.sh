#!/usr/bin/bash
#SBATCH --chdir=/home/grahman/projects/birdman-analyses-final
#SBATCH --output=/home/grahman/projects/birdman-analyses-final/slurm_out/obesity/%x.log
#SBATCH --partition=long
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -e
date; pwd; hostname

source ~/miniconda3/bin/activate birdman-analyses-final

CTX="Pick_closed-reference_OTUs-Greengenes-Illumina-16S-V4-5c6506"
QIITA_STR="qiita_study_id == 11548"
DIET_STR="diet in ('Regular chow', 'HFD')"
QUERY="where ${QIITA_STR} and ${DIET_STR}"

SAMP_FILE="data/obesity/intermediate/study_11548_samples.txt"
MD_FILE="data/obesity/intermediate/study_11548.tsv"
TBL_FILE="data/obesity/intermediate/study_11548.biom"

redbiom search metadata "${QUERY}" > $SAMP_FILE
redbiom fetch samples --from $SAMP_FILE --context $CTX --output $TBL_FILE
redbiom fetch sample-metadata --from $SAMP_FILE --context $CTX --output $MD_FILE
