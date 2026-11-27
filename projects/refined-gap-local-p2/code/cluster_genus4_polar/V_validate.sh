#!/bin/bash
#SBATCH --comment=validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=V_validate.%j.out
# V_validate.sh -- run every stage at genus 2 and 3 with 6 workers and compare with the stored results (ref/).
# All collectors print VALIDATION OK / VALIDATION FAIL lines; validate_summary.log must contain no FAIL before genus 4 is submitted.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/env.sh"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
NW=6
workers() { # workers <script.m> <prefix args...> [-- suffix]
  local script=$1; shift; local prefix=() suffix=() seen=0
  for a in "$@"; do if [ "$a" == "--" ]; then seen=1; continue; fi; if [ $seen -eq 0 ]; then prefix+=("$a"); else suffix+=("$a"); fi; done
  for w in $(seq 1 $NW); do wolframscript -file "$script" "${prefix[@]}" $w $NW "${suffix[@]}" > "worker_${script%.m}_w${w}.out" 2>&1 & done; wait; }
rm -f validate_summary.log
for g in 2 3; do
  K=$((3*g-3))
  echo "=== genus $g: local P2 ===" | tee -a validate_summary.log
  wolframscript -file lp2_prep.m $g
  workers lp2_sample.m $g $K
  wolframscript -file lp2_reshape.m $g $NW
  workers lp2_recon.m $g
  wolframscript -file lp2_collect.m $g
  echo "=== genus $g: resolved conifold ===" | tee -a validate_summary.log
  wolframscript -file rc_prep.m $g
  wolframscript -file rc_prep_x.m $g $K
  wolframscript -file rc_const.m $g $K
  workers rc_sample.m $g $K -- closed
  wolframscript -file rc_collect.m $g $K
  workers rc_polar_x.m $g
  wolframscript -file rc_polar_collect.m $g
  workers rc_polarCY_sample.m $g $K
  wolframscript -file rc_polarCY_reshape.m $g $NW
  workers rc_polarCY_fit.m $g
  wolframscript -file rc_polarCY_collect.m $g
  wolframscript -file compare_g4.m $g
done
wolframscript -file lp2_checks.m 3
wolframscript -file rc_check_ffc.m 9
grep -h "VALIDATION" lp2_collect.log rc_const.log rc_collect.log rc_polar_collect.log rc_polarCY_collect.log rc_check_ffc.log >> validate_summary.log
grep -h "anomaly\|GV d = 1..7\|gap and leading" lp2_checks_g3.log >> validate_summary.log
echo "=== summary ==="; cat validate_summary.log
if grep -q "VALIDATION FAIL" validate_summary.log; then echo "VALIDATION FAILED -- do not submit genus 4"; else echo "ALL VALIDATIONS OK -- ready for submit_genus4.sh"; fi
