#!/bin/bash
#SBATCH --comment=prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=S1_prep.%j.out
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
wolframscript -file lp2_prep.m 4
wolframscript -file rc_prep.m 4
wolframscript -file rc_prep_x.m 4 9
wolframscript -file rc_const.m 4 9
