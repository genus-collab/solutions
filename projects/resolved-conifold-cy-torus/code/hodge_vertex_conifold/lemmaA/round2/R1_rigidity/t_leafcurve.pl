#!/usr/bin/perl
# t_leafcurve.pl -- checks of the conifold leaf curve used for the closed-form (polar-order) relations.
# (1) A_j(t) from Theorem 11.1 of leaf_structure_proof.md agrees with the explicit q^(1), q^(2), q^(3) of its s.11.4.
# (2) After eliminating the lam_1 leaves by the string equation, the effective leaf coefficient is
#        b_j(t) = A_j(t) (-c^2/A_1(t))^j / c^2 = (-1)^j (Pi_j(f)^2/j) t^(j-1) rho_j(t),  rho_j = H_j(T)/H_1(T)^j,
#     H_j = 2F1(j x/c, j x/c; j; T), T = 1 - e^-t, framing f = -x/c;  and [t^1] log rho_j = 0, [t^2] log rho_j = -j(j-1) sigma^2/(4(j+1)).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my $N = 7; my $bad = 0; my $n = 0;
for my $x (R('2/7'), R('-3/5'), R('5/3'), R(3)) {
  my $c = 1 + $x; my $f = -$x/$c; my $sig = $f*(1+$f);
  my $Tt = T_of_t_series($N+1);
  # s* = 2/(c (1 + sqrt(1-u))), u = 4 x T/c^2 (series in T)
  my @u = (R(0), 4*$x/($c*$c), map { R(0) } 2..$N+1);
  my @omu = map { -$_ } @u; $omu[0] = R(1);
  my $sq = series_pow(\@omu, R('1/2'), $N+1);
  my @den = map { $sq->[$_] } 0..$N+1; $den[0] = $den[0] + 1; @den = map { $_ * $c } @den;
  my $sstar = series_inv(\@den, $N+1); $sstar = [ map { 2*$_ } @$sstar ];
  # dphi_T = s*/(1 - T s*)
  my @Ts = (R(0), @$sstar[0..$N]); my @omTs = map { -$_ } @Ts; $omTs[0] = R(1);
  my $dd = series_mul($sstar, series_inv(\@omTs, $N+1), $N+1);
  my @dphi = (R(0)); for my $k (1..$N+1) { $dphi[$k] = $dd->[$k-1]/$k; }
  my $e = series_exp([ map { $_ / $c } @dphi ], $N+1);
  my @zetaT = (R(0), @$e[0..$N]);            # zeta = T e^{dphi/c} as series in T
  my $zeta = series_compose(\@zetaT, $Tt, $N+1);   # in t
  my %A;
  for my $j (1..5) {
    my $Cj = $c*$c; my $q = R(1); for my $i (0..$j-2) { $q = $q * (1 - R($j)/$c + $i); } $Cj = $Cj * $q*$q / (fact($j)*fact($j-1));
    my $h = series_compose(hyp2f1_series($j*$x/$c, $j*$x/$c, R($j), $N+1), $Tt, $N+1);
    my $zj = R(1) ; my $zp = [ R(1), map { R(0) } 1..$N+1 ]; for (1..$j) { $zp = series_mul($zp, $zeta, $N+1); }
    my $prod = series_mul($zp, $h, $N+1);   # zeta^j 2F1, divisible by t
    $A{$j} = [ map { R((-1)**$j) * $Cj * $prod->[$_+1] } 0..$N ];   # A_j(t) = (-1)^j C_j t^-1 ...
  }
  # explicit values (s.11.4)
  my %ex = (
    '1,0' => -$c*$c, '1,1' => -(1-$x)**2/2, '2,1' => (1-$x)**2/2,
    '1,2' => -(2 - 4*$x + 3*$x**2 - 4*$x**3 + 2*$x**4)/(12*$c**2), '2,2' => 6*(1-$x)**4/(12*$c**2), '3,2' => -((2-$x)**2)*((1-2*$x)**2)/(12*$c**2),
    '1,3' => (1-$x)**2 * (-(3 + 5*$x**2 + 3*$x**4)/72)/$c**4, '2,3' => (1-$x)**2 * ((7 - 20*$x + 22*$x**2 - 20*$x**3 + 7*$x**4)/24)/$c**4,
    '3,3' => (1-$x)**2 * (-((2-$x)**2)*((1-2*$x)**2)/8)/$c**4, '4,3' => (1-$x)**2 * (((3-$x)**2)*((1-3*$x)**2)/36)/$c**4 );
  for my $k (sort keys %ex) { my ($j, $o) = split /,/, $k; $n++; if ($A{$j}[$o] != $ex{$k}) { $bad++; print "MISMATCH x=$x a_$j^($o): $A{$j}[$o] vs $ex{$k}\n"; } }
  # sum rule Sum_j A_j = -c^2 (orders <= 3 involve j <= 4)
  for my $o (0..3) { my $s = R(0); $s = $s + $A{$_}[$o] for 1..5; $n++; my $want = $o == 0 ? -$c*$c : R(0); if ($s != $want) { $bad++; print "sum rule fails at t^$o: $s\n"; } }
  # effective leaf coefficients
  my @alpha = map { -$_/($c*$c) } @{$A{1}};
  my $ia = series_inv(\@alpha, $N);
  my $h1 = series_compose(hyp2f1_series($x/$c, $x/$c, R(1), $N+1), $Tt, $N+1);
  for my $j (2..5) {
    my $iaj = [ R(1), map { R(0) } 1..$N ]; $iaj = series_mul($iaj, $ia, $N) for 1..$j;
    my $b = series_mul($A{$j}, $iaj, $N); $b = [ map { $_/($c*$c) } @$b ];
    my $hj = series_compose(hyp2f1_series($j*$x/$c, $j*$x/$c, R($j), $N+1), $Tt, $N+1);
    my $h1j = [ R(1), map { R(0) } 1..$N+1 ]; $h1j = series_mul($h1j, $h1, $N+1) for 1..$j;
    my $rho = series_mul($hj, series_inv($h1j, $N+1), $N+1);
    my $Pj = Pi_at($j, $f); my $top = R((-1)**$j) * $Pj*$Pj / $j;
    for my $o (0..$N-$j+1) { my $want = $o >= 0 ? $top * $rho->[$o] : R(0); $n++; if ($b->[$o + $j - 1] != $want) { $bad++; print "b_$j mismatch at t^", $o+$j-1, ": $b->[$o+$j-1] vs $want\n"; } }
    for my $o (0..$j-2) { $n++; unless ($b->[$o]->is_zero) { $bad++; print "b_$j has t^$o term\n"; } }
    my $lr = series_log($rho, $N);
    $n++; unless ($lr->[1]->is_zero) { $bad++; print "log rho_$j has a t^1 term $lr->[1]\n"; }
    my $w2 = R(-$j*($j-1)) * $sig*$sig / (4*($j+1)); $n++; if ($lr->[2] != $w2) { $bad++; print "[t^2] log rho_$j = $lr->[2], predicted $w2\n"; }
    if ($x == R('2/7')) { print "x=$x j=$j: log rho_j = ", join(', ', map { "$lr->[$_]" } 0..5), "\n"; }
  }
}
print "t_leafcurve: $n comparisons, $bad mismatches\n";
