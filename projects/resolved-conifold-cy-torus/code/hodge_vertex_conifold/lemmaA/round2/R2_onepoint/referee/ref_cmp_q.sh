#!/bin/bash
# Compare ref_xy_mod.pl (no engine table) with the direct evaluation of the Mathematica q-tables (ref_qtable.pl), mod p.
cd "$(dirname "$0")"
OUT=ref_cmp_q.out; : > $OUT
for g in 2 3 4; do for r in 1 2; do for t in 1/2 3 -2/7 5/3 4; do
  if [ $r -eq 1 ]; then J=$((3*g)); else J=5; fi
  perl ref_qtable.pl $t $g $r $J > tmpq_a.txt 2>&1
  perl ref_xy_mod.pl $t $g $r $J - > tmpq_b.txt 2>&1
  tot=0; bad=0; skip=0
  while read -r tup val; do
    line=$(grep -F "mu=($tup):" tmpq_b.txt | head -1)
    if echo "$line" | grep -q "Pi = 0"; then skip=$((skip+1)); continue; fi
    ref=$(echo "$line" | sed -E 's/.*: ref ([0-9]+) .*/\1/')
    tot=$((tot+1)); if [ "$ref" != "$val" ]; then bad=$((bad+1)); echo "  MISMATCH g=$g r=$r tau=$t mu=($tup): ref $ref data $val" >> $OUT; fi
  done < tmpq_a.txt
  echo "g=$g r=$r tau=$t: compared $tot, mismatches $bad, skipped (Pi=0) $skip" >> $OUT
done; done; done
echo "done" >> $OUT
