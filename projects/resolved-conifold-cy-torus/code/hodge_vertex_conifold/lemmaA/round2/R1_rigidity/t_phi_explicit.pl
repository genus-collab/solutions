#!/usr/bin/perl
# t_phi_explicit.pl -- the explicit flow phi_t(lambda) = alpha lambda - [Pcal_t]_+(lambda) + [Pcal_t]_+(0),
#   Pcal_t(lambda) = sum_{j>=2} Phat_j(t) u(lambda)^-j,  Phat_j = P_j alpha^j,  [u^-j]_+ = sum_{k=0}^{j} binom(j tau+k-1,k) lambda^(j-k).
# Tests (FC): log u_t(phi_t(lambda)) - log u(lambda) = -log C(t) (independent of lambda), i.e.  K_t(phi_t) - D = log C,
#   D(lambda) = tau log((phi-1)/(lambda-1)) - (1+tau) log(phi/lambda);  prints log C and compares with
#   log C_a = log(alpha e^{tau^2 t}/H_1(T)) (the critical-value form of C is checked in t_crit.pl).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M) = @ARGV; my $tau = R($taus); my $sig = $tau*(1+$tau); my $N = $M + 2;
sub gbin { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
my $Tt = T_of_t_series($N); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $N), $Tt, $N);
my %Pj; for my $j (2..$N+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $N), $Tt, $N); my $p = [R(1), map { R(0) } 1..$N]; $p = series_mul($p, $h1, $N) for 1..$j; my $rho = series_mul($hj, series_inv($p, $N), $N);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$N ]; }
my @al = (R(1), map { R(0) } 1..$N);
for my $it (0..$N+1) { my @nw = (R(1), map { R(0) } 1..$N); my $wp = [@al]; for my $j (2..$N+1) { $wp = series_mul($wp, \@al, $N); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp, $N); $nw[$_] = $nw[$_] + $q->[$_] for 0..$N; } @al = @nw; }
my %Ph; for my $j (2..$N+1) { my $ap = [R(1), map { R(0) } 1..$N]; $ap = series_mul($ap, \@al, $N) for 1..$j; $Ph{$j} = series_mul($Pj{$j}, $ap, $N); }
sub padd { my ($a, $b, $s) = @_; $s //= R(1); my $n = @$a > @$b ? @$a : @$b; return [ map { ($a->[$_] // R(0)) + $s*($b->[$_] // R(0)) } 0..$n-1 ]; }
sub pmul { my ($a, $b) = @_; my @c = map { R(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i]->is_zero; for my $j (0..$#$b) { $c[$i+$j] = $c[$i+$j] + $a->[$i]*$b->[$j]; } } return \@c; }
sub peval { my ($p, $x) = @_; my $v = R(0); my $pw = R(1); for my $c (@$p) { $v = $v + $c*$pw; $pw = $pw*$x; } return $v; }
# phi_t as t-series of polynomials
my @phi; for my $m (0..$M) { my @p = (R(0), $al[$m]); for my $j (2..$N+1) { my $c = $Ph{$j}[$m]; next if $c->is_zero; for my $k (0..$j-1) { $p[$j-$k] = ($p[$j-$k] // R(0)) - $c * gbin($j*$tau + $k - 1, $k); } } $phi[$m] = [ map { $_ // R(0) } @p ]; }
for my $m (1..$M) { print "xi_$m = ", join(' + ', map { "($phi[$m][$_])l^$_" } 0..$#{$phi[$m]}), "\n"; }
# check phi(1) = 1
for my $m (0..$M) { my $v = peval($phi[$m], R(1)); print "phi(1) at t^$m = $v\n" if $v != ($m == 0 ? 1 : 0); }
# t-series of polynomials ops
sub tmul { my ($A, $B) = @_; my @C = map { [R(0)] } 0..$M; for my $i (0..$M) { for my $j (0..$M-$i) { $C[$i+$j] = padd($C[$i+$j], pmul($A->[$i], $B->[$j])); } } return \@C; }
sub tlog1 { my ($X) = @_; # log(1 + X), X = O(t)
  my @res = map { [R(0)] } 0..$M; my $pw = [ [R(1)], map { [R(0)] } 1..$M ]; for my $k (1..$M) { $pw = tmul($pw, $X); for my $m (0..$M) { $res[$m] = padd($res[$m], $pw->[$m], R((-1)**($k+1))/$k); } } return \@res; }
# phi/lambda - 1 and (phi-1)/(lambda-1) - 1 as polynomials
my @X1 = ([R(0)]); my @X2 = ([R(0)]);
for my $m (1..$M) { my @p = @{$phi[$m]}; die "phi(0) != 0" unless $p[0]->is_zero; $X1[$m] = [ @p[1..$#p] ];
  # divide p by (lambda - 1): synthetic division
  my @q; my $r = R(0); for (my $i = $#p; $i >= 1; $i--) { $q[$i-1] = $p[$i] + $r; $r = $q[$i-1]; } die "phi(1) != 1 at $m" unless ($p[0] + $r)->is_zero; $X2[$m] = [ map { $_ // R(0) } @q ]; }
my $L1 = tlog1(\@X1); my $L2 = tlog1(\@X2);
my @Dt = map { padd([ map { $tau*$_ } @{$L2->[$_]} ], $L1->[$_], -(1+$tau)) } 0..$M;
# K_t(phi): sum_k kappa_k t^k [(1+tau) phi^k - tau (phi-1)^k]
my @Bn = (R(1)); for my $mm (1..$M+2) { my $s = R(0); for my $k (0..$mm-1) { $s = $s + binom($mm+1, $k) * $Bn[$k]; } $Bn[$mm] = -$s/($mm+1); }
my @kap = (R(0), $tau + R('1/2'), map { $Bn[$_]/($_*fact($_)) } 2..$M+2);
my @phit = @phi; my @phim1 = map { [ @$_ ] } @phi; $phim1[0] = [R(-1), R(1)];
my @KP = map { [R(0)] } 0..$M; my $pp = [ [R(1)], map { [R(0)] } 1..$M ]; my $pq = [ [R(1)], map { [R(0)] } 1..$M ];
for my $k (1..$M) { $pp = tmul($pp, \@phit); $pq = tmul($pq, \@phim1); for my $i (0..$M-$k) { $KP[$i+$k] = padd($KP[$i+$k], padd([ map { (1+$tau)*$_ } @{$pp->[$i]} ], $pq->[$i], -$tau), $kap[$k]); } }
my @logC; my $nonconst = 0;
for my $m (0..$M) { my $d = padd($KP[$m], $Dt[$m], R(-1)); my @nz = grep { !$d->[$_]->is_zero } 1..$#$d; $nonconst++ if @nz; push @logC, $d->[0] // R(0);
  print "t^$m: K_t(phi) - D = ", join(' + ', map { "($d->[$_])l^$_" } grep { !$d->[$_]->is_zero } 0..$#$d), "\n"; }
print "(FC) ", ($nonconst ? "FAILS ($nonconst orders non-constant)" : "holds: K_t(phi_t) - D is constant to order $M"), "\n";
# log C* and log C_a
my @et = map { R(1)/fact($_) } 0..$N+1; my @em1t = map { $et[$_+1] } 0..$N; my @lc = map { (1+$tau)*$et[$_] } 0..$N; $lc[0] = $lc[0] - $tau;
my $lcs = series_log(\@em1t, $M); my $llc = series_log(\@lc, $M); my @lCs = map { $lcs->[$_] + $tau*$llc->[$_] } 0..$M; $lCs[1] = $lCs[1] + $sig;
my $zt = series_mul(\@al, series_inv($h1, $N), $N); my $lz = series_log($zt, $M); my @lCa = map { $lz->[$_] } 0..$M; $lCa[1] = $lCa[1] + $tau*$tau;
print "log C     : ", join(' ', map { "$_" } @logC), "\n";
print "log C_a   : ", join(' ', map { "$_" } @lCa), "\n";
