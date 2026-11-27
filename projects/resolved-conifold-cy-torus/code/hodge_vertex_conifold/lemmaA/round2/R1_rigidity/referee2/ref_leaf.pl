#!/usr/bin/perl
# ref_leaf.pl -- referee 2: the leaf curve used in the polar relations, rebuilt mod p from leaf_structure_proof.md
# Theorem 11.1 / Corollary 11.2 (with its own zeta), independently of the author's r1lib.pm / t_leafcurve.pl.
#   x, c = 1+x, framing tau = -x/c, T = 1 - e^-t, u = 4xT/c^2, s* = 2/(c(1+sqrt(1-u))), U_1 = 1 - T s*,
#   zeta = c T s* U_1^(-(1-x)/c),  A_j = (-1)^j C_j t^-1 zeta^j 2F1(jx/c, jx/c; j; T),  C_j = c^2 [(1-j/c)_{j-1}]^2/(j!(j-1)!).
# Checks: (i) the explicit q^(1), q^(2), q^(3) of s.11.4; (ii) sum rule sum_j A_j = -c^2; (iii) b_j := A_j alpha^-j/c^2
# (alpha = -A_1/c^2) equals (-1)^j (Pi_j^2/j) t^(j-1) H_j/H_1^j; (iv) W := solution of W = 1 + sum_{j>=2} b_j W^j equals alpha.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rk;
my $N = 9; my ($n, $bad) = (0, 0);
for my $xs ('2/7', '-3/5', '5/3', '3', '-2/9') { my ($xn, $xd) = split m{/}, $xs; $xd //= 1; my $x = mfrac($xn, $xd); my $c = md(1 + $x); my $ic = mi($c);
  my $tau = md(-mm($x, $ic)); my $Tt = T_series($N + 1);
  # sqrt(1 - u) as a series in T: binomial series
  my @u = (0, mm(4, mm($x, mm($ic, $ic)))); my @sq = (1); my $g = 1; for my $k (1..$N+1) { $g = mm($g, mm(md(mi(2) - $k + 1), mi($k))); $sq[$k] = mm($g, mpow(md(-$u[1]), $k)); }
  my @den = map { mm($c, $_) } @sq; $den[0] = md($den[0] + $c); my $sst = ser_inv(\@den, $N+1); $sst = [ map { mm(2, $_) } @$sst ];   # s*(T)
  my @Ts = (0, @$sst[0..$N]); my @U1 = map { md(-$_) } @Ts; $U1[0] = 1;                                            # U_1 = 1 - T s*
  my $lU = ser_log(\@U1, $N+1); my $e = md(-mm(md(1 - $x), $ic)); my $pw = ser_exp([ map { mm($e, $_) } @$lU ], $N+1);
  my $zT = ser_mul(\@Ts, $pw, $N+1); $zT = [ map { mm($c, $_) } @$zT ];                                             # zeta as series in T
  my $zeta = ser_compose($zT, $Tt, $N+1);                                                                           # in t
  my %A; for my $j (1..6) { my $q = 1; for my $i (0..$j-2) { $q = mm($q, md(1 - mm($j, $ic) + $i)); } my $Cj = mm(mm($c, $c), mm(mm($q, $q), mi(mm(fct($j), fct($j-1)))));
    my $h = ser_compose(hyp_series(mm($j, mm($x, $ic)), mm($j, mm($x, $ic)), $j, $N+1), $Tt, $N+1);
    my $zp = [1, (0) x ($N+1)]; $zp = ser_mul($zp, $zeta, $N+1) for 1..$j; my $pr = ser_mul($zp, $h, $N+1);
    $A{$j} = [ map { mm(($j % 2) ? $P - 1 : 1, mm($Cj, $pr->[$_+1])) } 0..$N ]; }
  # explicit values s.11.4 (as rational functions of x)
  my $f = sub { my @co = @_; my $v = 0; my $p = 1; for my $cc (@co) { $v = md($v + mm($cc, $p)); $p = mm($p, $x); } return $v; };   # polynomial in x, ascending
  my $omx = md(1 - $x); my $ic2 = mm($ic, $ic); my $ic4 = mm($ic2, $ic2);
  my %ex = ('1,0' => md(-mm($c, $c)), '1,1' => md(-mm(mm($omx, $omx), mi(2))), '2,1' => mm(mm($omx, $omx), mi(2)),
    '1,2' => md(-mm($f->(2, -4, 3, -4, 2), mm(mi(12), $ic2))), '2,2' => mm(mm(6, mpow($omx, 4)), mm(mi(12), $ic2)),
    '3,2' => md(-mm(mm(mpow(md(2 - $x), 2), mpow(md(1 - 2*$x), 2)), mm(mi(12), $ic2))),
    '1,3' => mm(mm(mm($omx, $omx), md(-mm($f->(3, 0, 5, 0, 3), mi(72)))), $ic4),
    '2,3' => mm(mm(mm($omx, $omx), mm($f->(7, -20, 22, -20, 7), mi(24))), $ic4),
    '3,3' => mm(mm(mm($omx, $omx), md(-mm(mm(mpow(md(2 - $x), 2), mpow(md(1 - 2*$x), 2)), mi(8)))), $ic4),
    '4,3' => mm(mm(mm($omx, $omx), mm(mm(mpow(md(3 - $x), 2), mpow(md(1 - 3*$x), 2)), mi(36))), $ic4));
  for my $k (sort keys %ex) { my ($j, $o) = split /,/, $k; $n++; if ($A{$j}[$o] != $ex{$k}) { $bad++; print "x=$xs a_$j^($o) mismatch\n"; } }
  for my $o (0..4) { my $s = 0; $s = md($s + $A{$_}[$o]) for 1..6; $n++; my $want = $o == 0 ? md(-mm($c, $c)) : 0; if ($s != $want) { $bad++; print "x=$xs sum rule fails at t^$o\n"; } }
  my $alpha = [ map { md(-mm($_, $ic2)) } @{$A{1}} ]; my $ia = ser_inv($alpha, $N);
  my @Pi = (1, 1); for my $j (2..8) { my $v = 1; for my $b (1..$j-1) { $v = mm($v, md($j*$tau + $b)); } $Pi[$j] = mm($v, ifct($j-1)); }
  my $Pj = leaf_series($tau, \@Pi, $N, 6);
  for my $j (2..6) { my $iaj = [1, (0) x $N]; $iaj = ser_mul($iaj, $ia, $N) for 1..$j; my $b = ser_mul($A{$j}, $iaj, $N); $b = [ map { mm($_, $ic2) } @$b ];
    for my $o (0..$N-$j+1) { $n++; my $want = mm($Pi[$j], $Pj->{$j}[$o]); if ($b->[$o] != $want) { $bad++; print "x=$xs b_$j mismatch at t^$o\n"; } } }
  # W from the genus-0 equation (truncate j <= 6 valid to t^4)
  my $W = [1, (0) x $N]; for my $it (1..8) { my $s = [1, (0) x $N]; my $wp = [1, (0) x $N]; for my $j (1..6) { $wp = ser_mul($wp, $W, $N); next if $j < 2; my $tj = ser_mul($Pj->{$j}, $wp, $N); $s->[$_] = md($s->[$_] + mm($Pi[$j], $tj->[$_])) for 0..$N; } $W = $s; }
  for my $o (0..4) { $n++; if ($W->[$o] != $alpha->[$o]) { $bad++; print "x=$xs W != alpha at t^$o\n"; } } }
print "ref_leaf: $n comparisons, $bad mismatches\n";
