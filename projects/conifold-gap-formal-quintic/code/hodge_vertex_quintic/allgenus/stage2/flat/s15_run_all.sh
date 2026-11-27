#!/bin/sh
# s15_run_all.sh -- the full production of this folder (bash, Git Bash on Windows): 14 flat runs + 8 hatted runs, at most 4 processes,
# then the comparisons and the display table.  Per-run logs in logs/, JSON tables in $S15_SCRATCH (default: scratch/ in this folder).
# Wall/CPU times: logs/s15_run_all_times.log (one line per job; the jobs are single-threaded, so wall time = CPU time).
D="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "${S15_SCRATCH:-$D/scratch}"
J="${S15_SCRATCH:-$D/scratch}/jobs.txt"
: > "$J"
for p in 89 61; do
  echo "s15_hat.py $p 8 44 nat" >> "$J"; echo "s15_hat.py $p 6 34 nat" >> "$J"; echo "s15_hat.py $p 6 34 M0" >> "$J"
  for v in - nodil noT nosign stored; do echo "s15_flat.py $p 3 3 $v" >> "$J"; done
  echo "s15_flat.py $p 2 1" >> "$J"
done
echo "s15_hat.py 89 8 44 nat 1/3,-2/7" >> "$J"; echo "s15_hat.py 89 6 34 M0 1/3,-2/7" >> "$J"
echo "s15_flat.py 89 3 3 - all 1/3,-2/7" >> "$J"; echo "s15_flat.py 89 3 3 nodil all 1/3,-2/7" >> "$J"
T="$D/logs/s15_run_all_times.log"; echo "s15_run_all started $(date '+%Y-%m-%dT%H:%M:%S')" > "$T"
xargs -P 4 -I{} sh -c 't0=$(date +%s); py -I -B "'"$D"'"/$(echo {} | cut -d" " -f1) $(echo {} | cut -d" " -f2-) > /dev/null 2>&1; rc=$?; echo "$(( $(date +%s) - t0 )) s  rc=$rc  {}" >> "'"$T"'"' < "$J"
echo "s15_run_all finished $(date '+%Y-%m-%dT%H:%M:%S')" >> "$T"
py -I -B "$D/s15_compare.py" 89; py -I -B "$D/s15_compare.py" 61; py -I -B "$D/s15_compare.py" 89 _ab1o3_m2o7; py -I -B "$D/s15_rational.py" > /dev/null
