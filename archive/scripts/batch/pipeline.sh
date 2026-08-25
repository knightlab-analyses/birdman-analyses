#!/bin/bash

cd /home/grahman/projects/birdman-analyses-final

pwd; hostname; date

set -e

jid1=$(sbatch scripts/batch/0.01-simulate_data.sh | grep -oP "\d*")
jid2=$(sbatch --dependency=afterok:$jid1 scripts/batch/1.01-regress.sh | grep -oP "\d*")
