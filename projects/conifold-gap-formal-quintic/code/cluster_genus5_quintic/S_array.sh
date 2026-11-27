#!/bin/bash
#SBATCH --comment=g5array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm/S_array.%A_%a.out
# One array task of a genus-5 run. Usage (normally through submit_genus5.sh):
#   sbatch --array=0-<N-1> S_array.sh <N> --prime 89 --variant nat [--ab 1/3,-2/7] [--maxout 2]
# The array index must be 0-based. Resumable: a resubmitted task skips the inner graphs already in parts/<tag>/task<I>of<N>.json.
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
source ./env.sh || exit 1
NT="$1"; shift
$PY -I -B genus5_run.py --genus 5 "$@" --task "${SLURM_ARRAY_TASK_ID:?not an array job}" --ntasks "$NT"
