#!/usr/bin/perl
# t_res_ut.pl -- formal residue Res_mu u_t(mu)^n dmu (expansion at mu = infinity, order by order in t) versus the closed form
#   (-1)^(n-1) t^(n-1) e^{-n tau^2 t} Pi_n H_n(T),  H_n = 2F1(-n tau, -n tau; n; T), T = 1 - e^{-t}.
# u_t(mu)^n = u(mu)^n exp(-n K_t(mu)),  K_t(mu) = sum_k kappa_k t^k ((1+tau) mu^k - tau (mu-1)^k),  u^n = mu^-n (1 - 1/mu)^(n tau).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M, $NMAX) = @ARGV; my $tau = R($taus);
my @Bn = (R(1)); for my $m (1..$M+2) { my $s = R(0); for my $k (0..$m-1) { $s = $s + binom($m+1, $k) * $Bn[$k]; } $Bn[$m] = -$s/($m+1); }
my @kap = (R(0), $tau + R('1/2'), map { $Bn[$_]/($_*fact($_)) } 2..$M+2);
sub gbin { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
# K_t as t-series of polynomials in mu: K[k] = array of coeffs
my @K = ([R(0)]); for my $k (1..$M) { my @p = map { R(0) } 0..$k; $p[$k] = $p[$k] + (1+$tau); for my $i (0..$k) { $p[$i] = $p[$i] - $tau * binom($k, $i) * R((-1)**($k-$i)); } $K[$k] = [ map { $kap[$k]*$_ } @p ]; }
sub padd { my ($a, $b, $s) = @_; $s //= R(1); my $n = @$a > @$b ? @$a : @$b; return [ map { ($a->[$_] // R(0)) + $s*($b->[$_] // R(0)) } 0..$n-1 ]; }
sub pmul { my ($a, $b) = @_; my @c = map { R(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i]->is_zero; for my $j (0..$#$b) { $c[$i+$j] = $c[$i+$j] + $a->[$i]*$b->[$j]; } } return \@c; }
# exp(-n K_t) as t-series of polynomials
sub texp { my ($A) = @_; my @E = ([R(1)]); for my $m (1..$M) { my $s = [R(0)]; for my $k (1..$m) { $s = padd($s, pmul($A->[$k], $E[$m-$k]), R($k)); } $E[$m] = [ map { $_ / $m } @$s ]; } return \@E; }
my $Tt = T_of_t_series($M);
my $bad = 0;
for my $n (1..$NMAX) {
  my @A = map { [ map { -$n * $_ } @$_ ] } @K; my $E = texp(\@A);
  # Res mu^p * mu^-n (1-1/mu)^(n tau) = [mu^{n-p-1}] (1 - 1/mu)^(n tau) -> coefficient of mu^{-(p - n + 1)}: (-1)^q binom(n tau, q), q = p - n + 1
  my @res; for my $m (0..$M) { my $s = R(0); my $p = $E->[$m]; for my $i (0..$#$p) { my $q = $i - $n + 1; next if $q < 0; $s = $s + $p->[$i] * R((-1)**$q) * gbin($n*$tau, $q); } push @res, $s; }
  # closed form
  my $hn = series_compose(hyp2f1_series(-$n*$tau, -$n*$tau, R($n), $M), $Tt, $M);
  my @ex = map { (-$n*$tau*$tau)**$_ / fact($_) } 0..$M; my $cf = series_mul($hn, \@ex, $M);
  my @cfs = map { R(0) } 0..$M; for my $m (0..$M) { my $src = $m - ($n - 1); next if $src < 0; $cfs[$m] = R((-1)**($n-1)) * Pi_at($n, $tau) * $cf->[$src]; }
  my $ok = 1; for my $m (0..$M) { $ok = 0 if $res[$m] != $cfs[$m]; } $bad++ unless $ok;
  print "n=$n: Res u_t^n = (", join(', ', map { "$_" } @res), ")  closed form (", join(', ', map { "$_" } @cfs), ")  ", ($ok ? "OK" : "DIFF"), "\n"; }
print "t_res_ut: ", ($bad ? "$bad failures" : "all agree"), "\n";
