#!/bin/bash
# Referee sweep 1b: the tau = -2/7 runs of sweep 1 that stopped at Pi_7(-2/7) = 0 before the Pi = 0 guard was added.
cd "$(dirname "$0")"
E=../../../E_formulas
OUT=sweep1b.out
: > $OUT
run() { local t=$1 g=$2 n=$3 j=$4 tab=$5 pidx=${6:-0};
  timeout 3000 perl ref_xy_mod.pl "$t" $g $n $j "$tab" $pidx > tmp_sweep1b.txt 2>&1
  grep -E "MISMATCH [1-9]|INCONSISTENT|division|line [0-9]" tmp_sweep1b.txt | head -5 >> $OUT
  echo "Pi=0 lines: $(grep -c 'Pi = 0' tmp_sweep1b.txt)" >> $OUT
  tail -1 tmp_sweep1b.txt >> $OUT; }
for g in 3 4 5 6; do run -2/7 $g 1 $((3*g)) $E/engine_h${g}_r1.txt; done
run -2/7 7 1 21 ../engine_h7_r1.txt
for g in 2 3 4; do run -2/7 $g 2 7 $E/engine_h${g}_r2.txt; done
for g in 5 6; do run -2/7 $g 2 8 $E/engine_h${g}_r2.txt; done
echo "sweep1b done" >> $OUT
