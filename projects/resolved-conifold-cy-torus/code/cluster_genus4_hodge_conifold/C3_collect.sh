#!/bin/bash
#SBATCH --comment=hvccollect
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=C3_collect.%j.out
# C3_collect.sh -- collect the 120 chunk files: total versus the closed form F_4(0) + P_4 Li_{-5}(Q)/(eps1 eps2), gap, per stable graph
#                  (hvc_results_g4_Bm.m, hvc_table_g4_Bm.txt), then the structural analysis and the brief summary
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
wolframscript -file hvc_run.wls collect 4 B m && wolframscript -file hvc_analysis.wls 4 Bm && wolframscript -file hvc_summary.wls 4 Bm brief
