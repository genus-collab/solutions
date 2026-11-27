#!/bin/bash
#SBATCH --comment=g5validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/slurm/V_validate.%j.out
# Validation and calibration (about 10 minutes with 16 cores): genus 2-4 against the stored project results, genus-5 set-up, timing of
# three genus-5 graphs. Result: last line of logs/validate.log must be "G5_VALIDATE ALL CHECKS PASSED".
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}" || exit 1
source ./env.sh || exit 1
$PY -I -B genus5_run.py --validate --nproc "${SLURM_CPUS_PER_TASK:-16}"
