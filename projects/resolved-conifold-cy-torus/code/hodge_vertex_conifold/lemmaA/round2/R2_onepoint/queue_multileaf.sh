#!/bin/bash
# Sequential queue of multi-leaf checks of the x-y swap n-point formula against Route E's Hodge tables (one process at a time).
cd "$(dirname "$0")"
run() { local t=$1 g=$2 n=$3 j=$4; local out="runs_np_g${g}_n${n}_t${t}.out";
  timeout 7200 perl xy_npoint.pl $t $g $n $j > "$out" 2>&1; echo "g=$g n=$n tau=$t j<=$j: $(tail -1 "$out")"; }
run 1 3 3 4
run 2 3 3 4
run 1 4 3 4
run 2 4 3 3
run 1 3 4 3
run 1 5 3 4
run 1 4 4 3
run 1 2 5 3
echo "queue done"
