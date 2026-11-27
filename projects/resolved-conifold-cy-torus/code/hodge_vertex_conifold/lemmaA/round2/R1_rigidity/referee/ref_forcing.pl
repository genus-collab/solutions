#!/usr/bin/perl
# ref_forcing.pl -- referee's independent test of the FORCING THEOREM (s.9), without using T1, Lemma 8.1 or the
# odd-functional-equation argument.  At tau0 = -k0/q (mod p), genus G, odd level D:
#   unknowns y_D..y_{2G-3} (y_L = 0 for L < D: the (OL) hypothesis),
#   rows (i)  relation_k = 0, k = 1..2G-3, computed DIRECTLY: relation_k = [t^k] Sum_T prod P_t Y(T)/|Aut T|
#             with Y(T) from refkp's propagation and the leaf curve P_t (Theorem 11.1 + string);  [input (I1)]
#        (ii) two-leaf regularity at tau0: Y(m,n)(tau0) = 0 for q | n, m >= 1, m + n <= 2G-1 (m = 1: n y_{n-1} = 0).
# Reports whether y_D is forced to vanish, and whether ALL y_D..y_{2G-3} are forced (the proof gives F = 0 mod y^(2G-1)).
# Usage: perl ref_forcing.pl Gmin Gmax [all]   (all = also tuples outside the hypotheses)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use refkp;
my ($Gmin, $Gmax, $all) = @ARGV; set_mode($ENV{MODE} // 'mod'); set_prime($ENV{PRIME}) if $ENV{PRIME};
sub rank_mod { my @M = map { [ @$_ ] } @{$_[0]}; my $nc = $_[1]; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { unless (IZ($M[$i][$c])) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = IV($M[$rk][$c]); $M[$rk] = [ map { ML($_, $inv) } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || IZ($M[$i][$c]); my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = SB($M[$i][$j], ML($f, $M[$rk][$j])); } } $rk++; } return $rk; }
my %cache;
my ($nin, $nin_forced, $nin_allforced, $nout, $nout_forced) = (0, 0, 0, 0, 0);
for my $G ($Gmin..$Gmax) { my $K = 2*$G - 3;
  for my $q (4..2*$G-2) { for my $k0 (1..int(($q-1)/2)) { next unless gcd($k0, $q) == 1; next if $ENV{ONLYQ} && $q != $ENV{ONLYQ}; next if $ENV{ONLYK} && $k0 != $ENV{ONLYK};
    my $tau = NG(FR($k0, $q));
    my $S = kp_propagate($tau, $K); my $Pi = $S->{Pi};
    # relation rows
    my $N = $K; my $Tt = Tser($N); my $h1 = ser_comp(hyp2f1(NG($tau), NG($tau), ONE(), $N), $Tt, $N); my %rho;
    for my $j (2..$K+1) { my $hj = ser_comp(hyp2f1(NG(ML(NN($j), $tau)), NG(ML(NN($j), $tau)), NN($j), $N), $Tt, $N);
      my $pw = [ ONE(), map { ZERO() } 1..$N ]; $pw = ser_mul($pw, $h1, $N) for 1..$j; $rho{$j} = ser_mul($hj, ser_inv($pw, $N), $N); }
    my @rel = map { vzero($K) } 0..$K;
    for my $Dl (1..$K) { for my $T (parts_ge2($Dl, $Dl + 1)) { my $top = ONE(); my $pr = [ ONE(), map { ZERO() } 1..$K ];
        for my $t (@$T) { $top = ML($top, ML(NN(($t % 2) ? -1 : 1), ML($Pi->($t), IV(NN($t))))); $pr = ser_mul($pr, $rho{$t}, $K - $Dl); }
        $top = ML($top, IV(aut_f(@$T))); my $v = $S->{Y}{mkey(@$T)};
        for my $k ($Dl..$K) { vaxpy($rel[$k], ML($top, $pr->[$k - $Dl]), $v); } } }
    # two-leaf rows
    my @two;
    for (my $n = $q; $n <= 2*$G - 2; $n += $q) { for my $m (1..2*$G-1-$n) {
        if ($m == 1) { push @two, vscale(NN($n), vunit($K, $n - 2)); next; }
        my $v = $S->{Y}{mkey($m, $n)}; die "no Y($m,$n)" unless $v; push @two, $v; } }
    for (my $D = 3; $D <= $K; $D += 2) {
      my $hyp = ($k0 <= $D + 1 && $q >= $D + 2 && 3*($q - $k0) < 2*$G) ? 1 : 0;
      next unless $hyp || $all;
      my @cols = ($D-1 .. $K-1); my $nc = @cols;
      my @rows = map { my $r = $_; [ map { $r->[$_] } @cols ] } (@rel[1..$K], @two);
      my $r0 = rank_mod(\@rows, $nc);
      my @eD = map { ZERO() } 1..$nc; $eD[0] = ONE();
      my $r1 = rank_mod([ @rows, \@eD ], $nc);
      my $forced = ($r1 == $r0) ? 1 : 0; my $allf = ($r0 == $nc) ? 1 : 0;
      if ($hyp) { $nin++; $nin_forced += $forced; $nin_allforced += $allf; print "NOT FORCED (inside): G=$G D=$D q=$q k0=$k0\n" unless $forced; }
      else { $nout++; $nout_forced += $forced; print "  outside, forced: G=$G D=$D q=$q k0=$k0\n" if $forced && $ENV{VERBOSE}; }
    } } }
  print "G=$G done: inside so far $nin (forced $nin_forced, all y forced $nin_allforced); outside $nout (forced $nout_forced)\n";
}
print "ref_forcing G=$Gmin..$Gmax: inside the hypotheses $nin tuples, y_D forced in $nin_forced, all one-leaf values forced in $nin_allforced; outside $nout, forced $nout_forced\n";
