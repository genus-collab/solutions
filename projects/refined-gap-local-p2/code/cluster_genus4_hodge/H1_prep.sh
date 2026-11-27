#!/bin/bash
#SBATCH --comment=hvprep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=H1_prep.%j.out
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
wolframscript -file hv_prep.m 4
