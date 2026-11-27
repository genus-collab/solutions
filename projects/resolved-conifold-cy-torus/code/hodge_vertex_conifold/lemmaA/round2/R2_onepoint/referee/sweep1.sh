#!/bin/bash
# Referee sweep 1: ref_xy_mod.pl (formal Theorem 1 + P2, mod p) against Route E's TR-free engine tables, at framings the author did
# not use (non-integer tau and tau = 3..8), and at (g, n) the author did not reach (n = 4 at g = 4, 5; n = 5, 6, 7).
cd "$(dirname "$0")"
E=../../../E_formulas
OUT=sweep1.out
: > $OUT
run() { local t=$1 g=$2 n=$3 j=$4 tab=$5 pidx=${6:-0};
  timeout 3000 perl ref_xy_mod.pl "$t" $g $n $j "$tab" $pidx > tmp_sweep1.txt 2>&1
  grep -E "MISMATCH|no data|error|die|line [0-9]" tmp_sweep1.txt | head -5 >> $OUT
  tail -1 tmp_sweep1.txt >> $OUT; }
# n = 1
for g in 2 3 4 5 6; do for t in 1/2 -2/7 3/2 3 5 7 -5/3 1/9 11/4; do run $t $g 1 $((3*g)) $E/engine_h${g}_r1.txt; done; done
for t in 1/2 -2/7 3/2 5 6 7 8 -5/3; do run $t 7 1 21 ../engine_h7_r1.txt; done
# n = 2
for t in 1/2 -2/7 3/2 4 5 -5/3; do run $t 1 2 6 $E/engine_h1_r2.txt; done
for g in 2 3 4; do for t in 1/2 -2/7 3/2 4 5 -5/3; do run $t $g 2 7 $E/engine_h${g}_r2.txt; done; done
for g in 5 6; do for t in 1/2 -2/7 3 4 5; do run $t $g 2 8 $E/engine_h${g}_r2.txt; done; done
# n = 3
for t in 1/2 3 -2/7 5/3 4; do run $t 1 3 4 $E/engine_h1_r3.txt; done
for g in 2 3 4; do for t in 1/2 3 -2/7 5/3 4; do run $t $g 3 5 $E/engine_h${g}_r3.txt; done; done
for g in 5 6; do for t in 1/2 3 -2/7; do run $t $g 3 4 $E/engine_h${g}_r3.txt; done; done
# n = 4
for g in 2 3; do for t in 2 1/2 3 -2/7; do run $t $g 4 4 $E/engine_h${g}_r4.txt; done; done
for g in 4 5; do for t in 2 1/2 3; do run $t $g 4 3 $E/engine_h${g}_r4.txt; done; done
# n = 5
for g in 4 5; do for t in 2 1/2 3; do run $t $g 5 3 $E/engine_h${g}_r5.txt; done; done
# n = 6, 7
for t in 2 1/2; do run $t 5 6 2 $E/engine_h5_r6.txt; done
run 2 5 7 2 $E/engine_h5_r7.txt
# second prime on a subset
run 3 6 1 18 $E/engine_h6_r1.txt 1
run -2/7 7 1 21 ../engine_h7_r1.txt 1
run 4 5 2 8 $E/engine_h5_r2.txt 1
run 1/2 4 3 5 $E/engine_h4_r3.txt 1
run 3 4 4 3 $E/engine_h4_r4.txt 1
run 1/2 5 5 3 $E/engine_h5_r5.txt 1
echo "sweep1 done" >> $OUT
