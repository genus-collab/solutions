#!/bin/bash
#SBATCH --comment=hvvalidate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=V_validate.%j.out
# V_validate.sh -- run every stage at genus 2 and 3 with 6 workers and compare with the stored local runs (ref/), then the
# Hodge-integral / scalar-R-action checks up to genus 4. validate_summary.log must contain no VALIDATION FAIL before submit_genus4.sh.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
NW=6
workers() { # workers <script.m> <prefix args...>
  local script=$1; shift
  for w in $(seq 1 $NW); do wolframscript -file "$script" "$@" $w $NW > "worker_${script%.m}_w${w}.out" 2>&1 & done; wait; }
rm -f validate_summary.log
for g in 2 3; do
  K=$((3*g-3))
  echo "=== genus $g: Hodge vertices ===" | tee -a validate_summary.log
  wolframscript -file hv_prep.m $g
  workers hv_sample.m $g $K
  wolframscript -file hv_reshape.m $g $NW
  workers hv_recon.m $g
  wolframscript -file hv_collect.m $g
done
echo "=== Hodge integrals and scalar R-action (genus 4 PPZ on 379 graphs) ===" | tee -a validate_summary.log
wolframscript -file hv_checks.m 4
grep -h "VALIDATION" hv_reshape.log hv_collect.log hv_checks.log >> validate_summary.log
echo "=== summary ==="; cat validate_summary.log
if grep -q "VALIDATION FAIL" validate_summary.log; then echo "VALIDATION FAILED -- do not submit genus 4"; else echo "ALL VALIDATIONS OK -- ready for submit_genus4.sh"; fi
