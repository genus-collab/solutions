#!/bin/bash
#SBATCH --comment=rcCYfit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --array=1-20
#SBATCH --output=S9b_rc_polarCY_fit.%A_%a.out
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
export NTASKS=20
bash run_workers.sh rc_polarCY_fit.m 4
