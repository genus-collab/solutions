#!/bin/bash
#SBATCH --comment=hvsample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=72G
#SBATCH --time=12:00:00
#SBATCH --output=H2r_resume.%j.out
# H2r_resume.sh -- complete the genus-4 sampling: 9 workers (w of 9) walk through all 680 points and compute only those whose
# hv_samples_g4_p<idx>.mx is missing (about 54 after the run of 2026-10-29, in which 9 of the 120 kernels never started).
# About 6 missing points per worker at 270 s each: roughly half an hour.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
echo "sample files before: $(ls hv_samples_g4_p*.mx 2>/dev/null | wc -l) of 680"
export NTASKS=1
bash run_workers.sh hv_sample.m 4 9
status=$?
echo "sample files after: $(ls hv_samples_g4_p*.mx 2>/dev/null | wc -l) of 680 (launcher status $status)"
exit $status
