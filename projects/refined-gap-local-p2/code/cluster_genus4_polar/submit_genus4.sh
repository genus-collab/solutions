#!/bin/bash
# submit_genus4.sh -- submit the genus-4 stages with their dependencies (run from the directory, after V_validate.sh passed).
set -e
s1=$(sbatch --parsable S1_prep.sh)
s2=$(sbatch --parsable --dependency=afterok:$s1 S2_lp2_sample.sh)
s3a=$(sbatch --parsable --dependency=afterok:$s2 S3a_lp2_reshape.sh)
s3b=$(sbatch --parsable --dependency=afterok:$s3a S3b_lp2_recon.sh)
s4=$(sbatch --parsable --dependency=afterok:$s3b S4_lp2_collect.sh)
s5=$(sbatch --parsable --dependency=afterok:$s1 S5_rc_sample.sh)
s6=$(sbatch --parsable --dependency=afterok:$s5 S6_rc_collect.sh)
s7a=$(sbatch --parsable --dependency=afterok:$s1 S7a_rc_polar_x.sh)
s7b=$(sbatch --parsable --dependency=afterok:$s7a:$s6 S7b_rc_polar_collect.sh)
s8=$(sbatch --parsable --dependency=afterok:$s1 S8_rc_polarCY_sample.sh)
s9a=$(sbatch --parsable --dependency=afterok:$s8 S9a_rc_polarCY_reshape.sh)
s9b=$(sbatch --parsable --dependency=afterok:$s9a S9b_rc_polarCY_fit.sh)
s9c=$(sbatch --parsable --dependency=afterok:$s9b S9c_rc_polarCY_collect.sh)
s10=$(sbatch --parsable --dependency=afterok:$s4:$s7b:$s9c S10_compare.sh)
echo "submitted: S1 $s1, S2 $s2, S3a $s3a, S3b $s3b, S4 $s4, S5 $s5, S6 $s6, S7a $s7a, S7b $s7b, S8 $s8, S9a $s9a, S9b $s9b, S9c $s9c, S10 $s10"
