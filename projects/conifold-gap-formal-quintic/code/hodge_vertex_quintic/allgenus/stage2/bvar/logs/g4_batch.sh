#!/bin/sh
# genus-4 batch: sequential runs, each with a pool of 3 worker processes
D="$(cd "$(dirname "$0")/.." && pwd)"
S="$D/scratch"; mkdir -p "$S"
B="$D/bv_run.py"
cd "$S"
py -I -B "$B" run S2 4 89 --maxout 2 --nproc 3 --out "$S" > "$S/g4_S2_canon_out2.stdout" 2>&1
py -I -B "$B" run S2 4 89 --maxout 2 --norm stored --nproc 3 --out "$S" > "$S/g4_S2_stored_out2.stdout" 2>&1
py -I -B "$B" run S1 4 89 --maxout 2 --nproc 3 --out "$S" > "$S/g4_S1_canon_out2.stdout" 2>&1
py -I -B "$B" run S3 4 89 --maxout 2 --nproc 3 --out "$S" > "$S/g4_S3_canon_out2.stdout" 2>&1
py -I -B "$B" run S5 4 89 --maxout 2 --nproc 3 --out "$S" > "$S/g4_S5_canon_out2.stdout" 2>&1
py -I -B "$B" run S2 4 89 --nproc 3 --out "$S" > "$S/g4_S2_canon_full.stdout" 2>&1
echo G4_BATCH_DONE
