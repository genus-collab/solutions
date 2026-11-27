#!/bin/bash
# submit_resume.sh -- cancel the stalled genus-4 chain of 2026-10-29 (array tasks kept alive by hung wolframscript processes,
# and the stages waiting on them) and resubmit from the sampling resume:
#   H2r_resume (9 workers, missing points only) -> H3a reshape -> H3b reconstruction -> H4 collect,
# all with afterok dependencies and the robust launcher (run_workers.sh). Run from this directory.
set -e
cd "$(dirname "$0")"
for n in H2_sample.sh H3a_reshape.sh H3b_recon.sh H4_collect.sh H2r_resume.sh; do scancel -u "$USER" -n "$n" 2>/dev/null || true; done
sleep 5
echo "jobs still queued or running: $(squeue -u "$USER" -h | wc -l)"
echo "sample files present: $(ls hv_samples_g4_p*.mx 2>/dev/null | wc -l) of 680"
h2=$(sbatch --parsable H2r_resume.sh)
h3a=$(sbatch --parsable --dependency=afterok:$h2 H3a_reshape.sh)
h3b=$(sbatch --parsable --dependency=afterok:$h3a H3b_recon.sh)
h4=$(sbatch --parsable --dependency=afterok:$h3b H4_collect.sh)
echo "submitted: H2r $h2, H3a $h3a, H3b $h3b, H4 $h4"
