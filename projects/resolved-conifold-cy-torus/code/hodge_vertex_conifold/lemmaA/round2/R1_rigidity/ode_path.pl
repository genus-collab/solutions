#!/usr/bin/perl
# ode_path.pl -- test of the FLOW THEOREM along an arbitrary path P(t) of the slice (not the leaf curve):
#   the state of the genus-G linearised KP system (ode_state.pl equations) satisfies  Uhat_t(lambda) = [ Psi(phi_t(lambda)) ]_-,
#   phi_t(lambda) = W lambda - [Pcal_t]_+(lambda) + [Pcal_t]_+(0),  Pcal_t(lambda) = sum_j P_j(t) W^j u(lambda)^-j,
#   W = W(P(t)) the genus-0 scalar (W = 1 + sum_j Pi_j P_j W^j),  Psi(lambda) = sum_n y_{n-1} u(lambda)^n.
# Exact (Math::BigRat).  Path: P_j(t) = sum_{m >= j-1} c_{j,m} t^m with random small integers c_{j,m}.
# Usage: perl ode_path.pl K tau seed
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($K, $ts, $seed) = @ARGV; my $tau = R($ts); my $NMAX = $K + 2; my $KM = 2*$NMAX + 2; srand($seed // 1);
sub gbin { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
sub Pi { Pi_at($_[0], $tau) }
sub smul { series_mul($_[0], $_[1], $K) }
my (%Pj, %dPj);
for my $j (2..$NMAX) { $Pj{$j} = [ map { $_ >= $j-1 ? R(int(rand(7)) - 3) : R(0) } 0..$K ]; $dPj{$j} = [ map { ($_+1) * ($Pj{$j}[$_+1] // R(0)) } 0..$K ]; }
my @al = (R(1), map { R(0) } 1..$K);
for my $it (0..$K+1) { my @nw = (R(1), map { R(0) } 1..$K); my $wp = [@al]; for my $j (2..$NMAX) { $wp = smul($wp, \@al); my $q = smul([ map { Pi($j) * $_ } @{$Pj{$j}} ], $wp); $nw[$_] = $nw[$_] + $q->[$_] for 0..$K; } @al = @nw; }
my @alp = ([R(1), map { R(0) } 1..$K]); $alp[$_] = smul($alp[$_-1], \@al) for 1..2*$NMAX+2;
my $ali = series_inv(\@al, $K); my @alm = ([R(1), map { R(0) } 1..$K]); $alm[$_] = smul($alm[$_-1], $ali) for 1..2*$NMAX+2;
my (%c, %d); for my $n (1..$KM) { for my $i (1..$n) { $c{"$n,$i"} = R($i)/R($n) * gbin($tau*$n + $n - $i - 1, $n - $i); } }
for my $k (1..$KM) { for my $n (1..$k) { $d{"$k,$n"} = R((-1)**($k-$n)) * gbin($tau*$n, $k - $n); } }
my @y = (R(0), map { R(int(rand(19)) - 9) } 1..$KM);
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
# state in the y-frame: Uhat_t(y) = sum_n U_n alpha^-n u^(y)^n ; coefficient of y^kk valid for kk <= K+1-m
my %st; for my $m (0..$K) { for my $kk (1..$K+1-$m) { my $s = R(0); for my $n (1..$kk) { next if $n > $NMAX; for my $mp (0..$m) { $s = $s + $alm[$n][$mp]*$U[$n][$m-$mp]*$d{"$kk,$n"}; } } $st{"$m,$kk"} = $s; } }
# explicit phi_t
my %Ph; for my $j (2..$NMAX) { my $ap = [R(1), map { R(0) } 1..$K]; $ap = smul($ap, \@al) for 1..$j; $Ph{$j} = smul($Pj{$j}, $ap); }
my @xi; for my $m (1..$K) { my %p = (1 => $al[$m]); for my $j (2..$NMAX) { my $cc = $Ph{$j}[$m]; next if $cc->is_zero; for my $k (0..$j-1) { $p{$j-$k} = ($p{$j-$k} // R(0)) - $cc * gbin($j*$tau + $k - 1, $k); } } $xi[$m] = \%p; }
if ($ENV{PERTURB}) { $xi[2]{2} = ($xi[2]{2} // R(0)) + R(1); }   # sensitivity test
# Psi(phi_t(lambda)) as Laurent series in lambda, via Taylor expansion (cf. fit_flow.pl)
my $NN = $K + 6; my $u = series_mul([R(0), R(1)], series_pow([R(1), R(-1)], $tau, $NN), $NN);
my $psi = [R(0), map { R(0) } 1..$NN]; my $up = [R(1), map { R(0) } 1..$NN]; for my $n (1..$NN) { $up = series_mul($up, $u, $NN); next if $n < 2 || $n > $KM; $psi->[$_] = $psi->[$_] + $y[$n-1]*$up->[$_] for 0..$NN; }
my %psiL; for my $n (0..$NN) { $psiL{-$n} = $psi->[$n] unless $psi->[$n]->is_zero; }
sub lmul { my ($a, $b, $lo) = @_; my %c; for my $i (keys %$a) { for my $j (keys %$b) { my $k = $i + $j; next if $k < $lo; $c{$k} = ($c{$k} // R(0)) + $a->{$i}*$b->{$j}; } } return \%c; }
sub ladd { my ($a, $b, $s) = @_; my %c = %$a; for my $k (keys %$b) { $c{$k} = ($c{$k} // R(0)) + ($s // R(1))*$b->{$k}; } return \%c; }
sub ld { my $a = shift; my %c; for my $k (keys %$a) { next if $k == 0; $c{$k-1} = ($c{$k-1} // R(0)) + R($k)*$a->{$k}; } return \%c; }
my $lo = -($K + 3);
my @res = map { {} } 0..$K; $res[0] = { %psiL };
my @dpow = ([ {0 => R(1)}, map { {} } 1..$K ]); my @der = (\%psiL); for my $k (1..$K) { $der[$k] = ld($der[$k-1]); }
for my $k (1..$K) { my @prev = @{$dpow[$k-1]}; my @cur = map { {} } 0..$K;
  for my $i (0..$K) { next unless %{$prev[$i]}; for my $m (1..$K-$i) { $cur[$i+$m] = ladd($cur[$i+$m], lmul($prev[$i], $xi[$m], $lo - 60)); } }
  $dpow[$k] = \@cur; my $f = R(1)/fact($k); for my $i (0..$K) { next unless %{$cur[$i]}; $res[$i] = ladd($res[$i], lmul($cur[$i], $der[$k], $lo), $f); } }
my ($n, $bad) = (0, 0);
for my $m (0..$K) { for my $kk (1..$K+1-$m) { $n++; my $a = $st{"$m,$kk"}; my $b = $res[$m]{-$kk} // R(0); $bad++ if $a != $b; } }
print "ode_path K=$K tau=$ts seed=$seed: $n state coefficients compared with [Psi(phi_t)]_-, $bad mismatches\n";
