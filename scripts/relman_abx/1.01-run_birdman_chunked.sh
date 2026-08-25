#!/bin/bash
#SBATCH --job-name=birdman_relman_abx
#SBATCH --output=slurm_out/relman_abx/%x.%a.out
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --array=1-20
#SBATCH --export=ALL

# fits one birdman model per feature, split across the job array.
#
# each feature is fit independently, so the array size only changes how the
# work is divided -- not the results. output filenames are indexed against the
# full table, so they are stable no matter what --array you use.
#
# submit from the repo root, after creating the log dir. pass the output dir
# explicitly rather than relying on the environment propagating:
#   mkdir -p slurm_out/relman_abx
#   sbatch --export=ALL,BIRDMAN_INFERENCE_DIR=/path/to/inferences \
#       scripts/relman_abx/1.01-run_birdman_chunked.sh
# defaults to <repo>/inferences/relman_abx if unset.

set -e
pwd; hostname; date

# under sbatch the script runs as a copy in /var/spool, so BASH_SOURCE is
# useless here -- SLURM_SUBMIT_DIR is the directory sbatch was called from
PROJ_DIR="${BIRDMAN_ANALYSES_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}}"

if [ ! -f "$PROJ_DIR/src/relman_abx/run_birdman_chunked.py" ]; then
    echo "ERROR: PROJ_DIR=$PROJ_DIR is not the birdman-analyses repo root." >&2
    echo "Submit from the repo root, or pass BIRDMAN_ANALYSES_DIR=/path/to/repo." >&2
    exit 1
fi
cd "$PROJ_DIR"

if [ -z "${CONDA_BASE:-}" ]; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
fi
if [ ! -f "${CONDA_BASE:-}/etc/profile.d/conda.sh" ]; then
    for c in "$HOME/software/miniconda3" "$HOME/miniconda3" "$HOME/anaconda3" /opt/conda; do
        if [ -f "$c/etc/profile.d/conda.sh" ]; then CONDA_BASE="$c"; break; fi
    done
fi
if [ ! -f "${CONDA_BASE:-}/etc/profile.d/conda.sh" ]; then
    echo "ERROR: no conda install found. set CONDA_BASE=/path/to/conda." >&2
    exit 1
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "${BIRDMAN_CONDA_ENV:-q2-birdman-dev-7}"

# cmdstan ships inside the conda env; PYTHONNOUSERSITE keeps a stray ~/.local
# install from shadowing the env's packages
export CMDSTAN="${CMDSTAN:-$CONDA_PREFIX/bin/cmdstan}"
export PYTHONNOUSERSITE=1
export PYTHONPATH="$PROJ_DIR:$PYTHONPATH"

# sbatch does not always propagate the submitting shell's environment, so this
# needs a working default rather than a hard requirement
OUTDIR="${BIRDMAN_INFERENCE_DIR:-$PROJ_DIR/inferences/relman_abx}"
LOGDIR="$PROJ_DIR/results/relman_abx/log"
mkdir -p "$OUTDIR" "$LOGDIR"

echo "chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX -> $OUTDIR"

time python src/relman_abx/run_birdman_chunked.py \
    --inference-dir "$OUTDIR" \
    --num-chunks "$SLURM_ARRAY_TASK_MAX" \
    --chunk-num "$SLURM_ARRAY_TASK_ID" \
    --chains 4 \
    --num-iter 500 \
    --num-warmup 1000 \
    --beta-prior 10.0 \
    --disp-scale 0.5 \
    --re-prior 3.0 \
    --logfile "$LOGDIR/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!
