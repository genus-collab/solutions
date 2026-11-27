#!/bin/bash
#SBATCH --comment=lp2sample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --array=1-20
#SBATCH --output=S2_lp2_sample.%A_%a.out
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
export NTASKS=20
bash run_workers.sh lp2_sample.m 4 9
