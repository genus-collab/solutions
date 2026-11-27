#!/bin/bash
#SBATCH --comment=lp2reshape
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=S3a_lp2_reshape.%j.out
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
wolframscript -file lp2_reshape.m 4 120
