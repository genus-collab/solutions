#!/bin/sh
# Referee 2 (route TR2): commands behind out/*.txt (run from this folder; at most 4 at a time).
# Engine sanity tests (seconds):
timeout 900 perl t_engine.pl
# Lemma D / Theorem 1 directly in TR (independent engine RefTR.pm), all (g,n) with n <= NM + G - g, two primes:
timeout 3000 perl lemmaD_ref.pl 3 1 8 1 2          > out/lemmaD_tau3_1_G8.txt
timeout 3000 perl lemmaD_ref.pl 5 2 8 1 2          > out/lemmaD_tau5_2_G8.txt
timeout 3000 perl lemmaD_ref.pl -3 4 8 1 2         > out/lemmaD_tau-3_4_G8.txt
timeout 3000 perl lemmaD_ref.pl 1 3 8 1 2          > out/lemmaD_tau1_3_G8.txt
timeout 3000 perl lemmaD_ref.pl 123456789 987654321 8 1 2 > out/lemmaD_tauRandom_G8.txt
timeout 3600 perl lemmaD_ref.pl 7 1 9 1 2          > out/lemmaD_tau7_1_G9.txt
timeout 3600 perl lemmaD_ref.pl -7 3 9 1 2         > out/lemmaD_tau-7_3_G9.txt
# closed formulas (5.2), (6.1), (5.4), Norlund, general-n (n = 2, 3, 4) against TR:
timeout 3600 perl closed_ref.pl 3 1 8 6 8 0        > out/closed_tau3_1_G8.txt
timeout 3600 perl closed_ref.pl -7 3 6 4 6 4       > out/closed_tau-7_3_G6_n3.txt
( timeout 3000 perl closed_ref.pl 123456789 987654321 7 6 8 0; timeout 3000 perl closed_n34.pl 123456789 987654321 4 4 2 3 ) > out/closed_tauRandom.txt
# P2 against Route E's Hodge engine tables (multi-leaf):
for hr in "4 4" "4 5" "5 2" "5 3" "6 2" "5 4" "6 3"; do timeout 3000 perl p2_engine_ref.pl 2 5 $hr 2; done > out/p2_engine_tau2_5.txt
for hr in "2 3" "3 4" "4 4" "4 5" "5 2" "5 3" "6 2" "5 4" "6 3"; do timeout 3000 perl p2_engine_ref.pl 123456789 987654321 $hr 2; done > out/p2_engine_tauRandom.txt
# claim 9 (tau = 1 one-line identity):
timeout 900 perl tau1_identity.pl 7 34 > out/tau1_identity.txt
