#!/bin/bash
#SBATCH --comment=hvcorbits
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=72G
#SBATCH --time=48:00:00
#SBATCH --array=1-20
#SBATCH --output=C2_orbits.%A_%a.out
# C2_orbits.sh -- the 27172 genus-4 decorated graphs (mode B) in 120 chunks (20 array tasks x 6 single-kernel workers), exact terms
#                 in Q and x and polar parts at Q = 1; worker w of 120 computes chunk w and appends every finished decorated graph to
#                 hvc_terms_g4_Bm_c<w>of120.m at once, so a resubmission only computes what is missing.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
export NTASKS=20
bash run_workers.sh hvc_run.wls orbits 4 B m
