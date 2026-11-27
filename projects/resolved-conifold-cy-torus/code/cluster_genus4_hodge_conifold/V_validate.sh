#!/bin/bash
#SBATCH --comment=hvcvalidate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=V_validate.%j.out
# V_validate.sh -- checks (Hodge integrals, PPZ scalar action, rank-4 triviality, genus-2 sign test), then every stage at genus 2 and 3
# with 6 workers (prep, orbits in 6 chunks, collect, analysis) compared with the stored local runs in ref/.
# validate_summary.log must contain no VALIDATION FAIL before submit_genus4.sh.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
NW=6
rm -f validate_summary.log hvc_validate.log
echo "=== checks ===" | tee -a validate_summary.log
wolframscript -file hvc_run.wls checks | grep -E "CHECKS|trivial:" | tee -a validate_summary.log
wolframscript -file hvc_rank4_check.wls | tail -1 | tee -a validate_summary.log
wolframscript -file hvc_run.wls g2sign | grep "sign test" | tee -a validate_summary.log
for g in 2 3; do
  echo "=== genus $g, mode B ===" | tee -a validate_summary.log
  wolframscript -file hvc_run.wls prep $g B m
  for w in $(seq 1 $NW); do wolframscript -file hvc_run.wls orbits $g B m $w $NW > "worker_validate_g${g}_w${w}.out" 2>&1 & done; wait
  wolframscript -file hvc_run.wls collect $g B m
  wolframscript -file hvc_analysis.wls $g Bm > /dev/null
  wolframscript -file hvc_validate.wls $g
done
grep -h "VALIDATION" hvc.log hvc_validate.log >> validate_summary.log
echo "=== summary ==="; cat validate_summary.log
if grep -q "VALIDATION FAIL" validate_summary.log; then echo "VALIDATION FAILED -- do not submit genus 4"; else echo "ALL VALIDATIONS OK -- ready for submit_genus4.sh"; fi
