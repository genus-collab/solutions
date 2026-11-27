#!/bin/bash
#SBATCH --comment=hvreshape
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=H3a_reshape.%j.out
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
wolframscript -file hv_reshape.m 4 120
