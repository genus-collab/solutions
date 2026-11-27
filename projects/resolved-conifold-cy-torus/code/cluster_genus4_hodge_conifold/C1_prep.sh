#!/bin/bash
#SBATCH --comment=hvcprep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=C1_prep.%j.out
# C1_prep.sh -- genus-4 factors of R' = R_Z h_w^-1 (mode B) with Q and x symbolic to z^9 -> hvc_factors_g4_Bm.mx,
#               and the decorated graphs with Hodge vertices -> hvc_orbits_g4_hodge.mx (one kernel)
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
wolframscript -file hvc_run.wls prep 4 B m
