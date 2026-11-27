#!/bin/bash
# submit_genus4.sh -- submit the genus-4 Hodge-vertex stages with their dependencies (run from this directory after V_validate.sh passed).
set -e
h1=$(sbatch --parsable H1_prep.sh)
h2=$(sbatch --parsable --dependency=afterok:$h1 H2_sample.sh)
h3a=$(sbatch --parsable --dependency=afterok:$h2 H3a_reshape.sh)
h3b=$(sbatch --parsable --dependency=afterok:$h3a H3b_recon.sh)
h4=$(sbatch --parsable --dependency=afterok:$h3b H4_collect.sh)
echo "submitted: H1 $h1, H2 $h2, H3a $h3a, H3b $h3b, H4 $h4"
