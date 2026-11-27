#!/bin/bash
#SBATCH --comment=g5collect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/C_collect.%j.out
# Merge the part files of one run and check it:  sbatch C_collect.sh --prime 89 --variant nat [--ab ..] [--maxout ..]
# Result: results/<tag>.json and logs/<tag>/collect.log (last line "GENUS5_RUN COLLECT (<tag>) ALL CHECKS PASSED" or the failures).
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
source ./env.sh || exit 1
$PY -I -B genus5_run.py --genus 5 "$@" --collect
