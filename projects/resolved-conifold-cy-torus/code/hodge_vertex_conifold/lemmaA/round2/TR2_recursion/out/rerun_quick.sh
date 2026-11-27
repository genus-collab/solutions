#!/bin/sh
# re-runs the quick checks of route TR2 and stores their outputs (run from TR2_recursion)
timeout 900 perl onepoint.pl 2 5 3 > out/onepoint_tau2_g5.txt 2>&1
timeout 900 perl onepoint.pl 3/7 4 3 > out/onepoint_tau3_7_g4.txt 2>&1
timeout 900 perl onepoint.pl -5/3 4 3 > out/onepoint_taum5_3_g4.txt 2>&1
timeout 900 perl closed_vs_61.pl -2/7 6 > out/closed_vs_61_taum2_7.txt 2>&1
( timeout 600 perl norlund_check.pl 3/7 6 6; timeout 600 perl norlund_check.pl -13/5 5 6; timeout 600 perl norlund_check.pl 2 5 6 ) > out/norlund_check.txt 2>&1
( timeout 600 perl trace_qle.pl 2 1 4; timeout 600 perl trace_qle.pl 3 7 3; timeout 600 perl trace_qle.pl 1 1 4 ) > out/trace_qle.txt 2>&1
( for g in 1 2 3; do timeout 600 perl xy_point.pl 2 1 $g 1 1; done
  for t in "-1 2" "-2 7" "7 1" "1 5"; do timeout 900 perl xy_point.pl $t 4 1 31; done
  for g in 1 2 3; do timeout 600 perl xy_point.pl 3 7 $g 2 5; done
  timeout 1500 perl xy_point.pl -5 3 4 2 11; timeout 1500 perl xy_point.pl 2 1 4 2 3; timeout 900 perl xy_point.pl 1 2 2 2 4; timeout 900 perl xy_point.pl -1 2 3 2 33
  for g in 0 1 2; do timeout 900 perl xy_point.pl -5 3 $g 3 7; done
  timeout 3000 perl xy_point.pl 3 7 3 3 13
  timeout 1200 perl xy_point.pl 2 1 0 4 21; timeout 2400 perl xy_point.pl 2 1 1 4 17; timeout 2400 perl xy_point.pl 3 7 1 4 22 ) > out/xy_point_all.txt 2>&1
echo finished > out/rerun_quick.done
