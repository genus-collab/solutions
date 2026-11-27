#!/usr/bin/perl
# ode_state.pl -- exact (Math::BigRat) version of ode_leaf.pl that prints the state along the leaf curve in the rescaled y-frame:
#   Uhat_t(y) := U(alpha/u^(y); t) = sum_n U_n(t) alpha^{-n} u^(y)^n,   u^(y) = y(1-y)^tau,
#   F_t(y)    := Vhat_t(u^(y))  (the Hankel datum of the two-point function at P(t)).
# Usage: perl ode_state.pl K tau L0      (basic solution: y_{L0} = 1, all other one-leaf values 0)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($K, $ts, $L0) = @ARGV; my $tau = R($ts); my $sig = $tau*(1+$tau); my $NMAX = $K + 2; my $KM = 2*$NMAX + 2;
sub gbin { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
sub Pi { Pi_at($_[0], $tau) }
sub smul { series_mul($_[0], $_[1], $K) }
my $Tt = T_of_t_series($K); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $K), $Tt, $K);
my (%Pj, %dPj);
for my $j (2..$NMAX) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $K), $Tt, $K); my $p = [R(1), map { R(0) } 1..$K]; $p = smul($p, $h1) for 1..$j; my $rho = smul($hj, series_inv($p, $K));
  my $cc = R((-1)**$j) * Pi($j) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$K ]; $dPj{$j} = [ map { ($_+1) * ($Pj{$j}[$_+1] // R(0)) } 0..$K ]; }
my @al = (R(1), map { R(0) } 1..$K);
for my $it (0..$K+1) { my @nw = (R(1), map { R(0) } 1..$K); my $wp = [@al]; for my $j (2..$NMAX) { $wp = smul($wp, \@al); my $q = smul([ map { Pi($j) * $_ } @{$Pj{$j}} ], $wp); $nw[$_] = $nw[$_] + $q->[$_] for 0..$K; } @al = @nw; }
my @alp = ([R(1), map { R(0) } 1..$K]); $alp[$_] = smul($alp[$_-1], \@al) for 1..2*$NMAX+2;
my $ali = series_inv(\@al, $K); my @alm = ([R(1), map { R(0) } 1..$K]); $alm[$_] = smul($alm[$_-1], $ali) for 1..2*$NMAX+2;
my (%c, %d); for my $n (1..$KM) { for my $i (1..$n) { $c{"$n,$i"} = R($i)/R($n) * gbin($tau*$n + $n - $i - 1, $n - $i); } }
for my $k (1..$KM) { for my $n (1..$k) { $d{"$k,$n"} = R((-1)**($k-$n)) * gbin($tau*$n, $k - $n); } }
my @y = map { R(0) } 0..$KM; $y[$L0] = R(1);
my (@U, @V, @F, %G, %GH); for my $n (1..$NMAX) { $U[$n] = [ map { R(0) } 0..$K ]; $V[$n] = [ map { R(0) } 0..$K ]; } for my $k (1..$KM) { $F[$k] = [ map { R(0) } 0..$K ]; }
for my $n (2..$NMAX) { $U[$n][0] = $y[$n-1]; }
for my $m (0..$K) {
  for my $n (1..$NMAX) { next if ($n-1)+$m > $K+1; my $v = R($n) * $U[$n][$m];
    for my $j (2..$NMAX) { for my $mp ($j-1..$m) { my $pj = $Pj{$j}[$mp]; next if $pj->is_zero; my $g = $G{"$j,$n"}; next unless $g && defined $g->[$m-$mp]; $v = $v + R($j)*$pj*$g->[$m-$mp]; } }
    $V[$n][$m] = $v; }
  for my $k (1..$KM) { next if ($k-1)+$m > 2*$K+2; my $s = R(0); for my $n (1..($k < $NMAX ? $k : $NMAX)) { my $dk = $d{"$k,$n"}; next if $dk->is_zero; for my $mp (0..$m) { next if $V[$n][$m-$mp]->is_zero; $s = $s + $alm[$n][$mp]*$dk*$V[$n][$m-$mp]; } } $F[$k][$m] = $s; }
  for my $a (1..$NMAX) { for my $b (1..$NMAX) { next if ($a-1)+($b-1)+$m > $K+1; my $gh = R(0);
      for my $j (1..$a) { my $ca = $c{"$a,$j"}; next if $ca->is_zero; for my $i (1..$b) { my $cb = $c{"$b,$i"}; next if $cb->is_zero; $gh = $gh + $ca*$cb*$F[$i+$j-1][$m]; } }
      $GH{"$a,$b"}[$m] = $gh; my $s = R(0); my $e = $a + $b - 1; for my $mp (0..$m) { next unless defined $GH{"$a,$b"}[$m-$mp]; $s = $s + $alp[$e][$mp]*$GH{"$a,$b"}[$m-$mp]; } $G{"$a,$b"}[$m] = $s; } }
  last if $m == $K;
  for my $n (1..$NMAX) { next if ($n-1)+$m+1 > $K+1; my $s = R(0);
    for my $i (2..$NMAX) { for my $mp ($i-2..$m) { my $dp = $dPj{$i}[$mp]; next if $dp->is_zero; my $g = $G{"$i,$n"}; next unless $g && defined $g->[$m-$mp]; $s = $s + $dp*$g->[$m-$mp]; } }
    $U[$n][$m+1] = $s / ($m+1); } }
# rescaled one-point function in y:  Uhat_t(y) = sum_n U_n alpha^{-n} u^(y)^n ;  F_t(y) = sum_k F_k y^k
for my $m (0..$K) { my @uh; for my $kk (1..$K+1+($ENV{EXTRA}//0)) { my $s = R(0); for my $n (1..$kk) { next if $n > $NMAX; for my $mp (0..$m) { $s = $s + $alm[$n][$mp]*$U[$n][$m-$mp]*$d{"$kk,$n"}; } } push @uh, $s; }
  print "t^$m Uhat(y): ", join('  ', map { "y$_:" . $uh[$_-1] } grep { !$uh[$_-1]->is_zero } 1..$K+1), "\n"; }
for my $m (0..$K) { print "t^$m F(y): ", join('  ', map { "y$_:" . $F[$_][$m] } grep { !$F[$_][$m]->is_zero } 1..$K+1), "\n"; }
