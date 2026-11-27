#!/bin/bash
# submit_genus4.sh -- submit the genus-4 stages with their dependencies (run from this directory after V_validate.sh passed).
set -e
cd "$(dirname "$0")"
c1=$(sbatch --parsable C1_prep.sh)
c2=$(sbatch --parsable --dependency=afterok:$c1 C2_orbits.sh)
c3=$(sbatch --parsable --dependency=afterok:$c2 C3_collect.sh)
echo "submitted: C1 $c1, C2 $c2, C3 $c3"
