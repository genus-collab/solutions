#!/bin/bash
#SBATCH --comment=rcsample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --array=1-20
#SBATCH --output=S5_rc_sample.%A_%a.out
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
export NTASKS=20
bash run_workers.sh rc_sample.m 4 9 -- closed
