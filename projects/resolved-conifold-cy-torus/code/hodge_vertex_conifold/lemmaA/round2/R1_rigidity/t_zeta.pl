#!/usr/bin/perl
# t_zeta.pl -- compare alpha(t) (genus-0 iteration on the leaf curve) with zeta H_1/t for two candidate zetas.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M) = @ARGV; my $tau = R($taus); my $N = $M + 2;
my $Tt = T_of_t_series($N); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $N), $Tt, $N);
my %Pj; for my $j (2..$N+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $N), $Tt, $N); my $p = [R(1), map { R(0) } 1..$N]; $p = series_mul($p, $h1, $N) for 1..$j; my $rho = series_mul($hj, series_inv($p, $N), $N);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$N ]; }
my @al = (R(1), map { R(0) } 1..$N);
for my $it (0..$N+1) { my @nw = (R(1), map { R(0) } 1..$N); my $wp = [@al]; for my $j (2..$N+1) { $wp = series_mul($wp, \@al, $N); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp, $N); $nw[$_] = $nw[$_] + $q->[$_] for 0..$N; } @al = @nw; }
print "alpha: ", join(' ', map { "$_" } @al[0..$M]), "\n";
# zeta/t from alpha: zeta/t = alpha/H_1
my $zt = series_mul(\@al, series_inv($h1, $N), $N); print "alpha/H1 = zeta/t: ", join(' ', map { "$_" } @$zt[0..$M]), "\n";
# candidate: zeta_alt/t = ((e^t-1)/t) e^{tau t} ((1+tau) e^t - tau)^tau
my @et = map { R(1)/fact($_) } 0..$N+1; my @em1t = map { $et[$_+1] } 0..$N;   # (e^t - 1)/t
my @etau = map { $tau**$_ / fact($_) } 0..$N;
my @lc = map { (1+$tau)*$et[$_] } 0..$N; $lc[0] = $lc[0] - $tau; my $lcp = series_pow(\@lc, $tau, $N);
my $alt = series_mul(series_mul(\@em1t, \@etau, $N), $lcp, $N); print "zeta_alt/t: ", join(' ', map { "$_" } @$alt[0..$M]), "\n";
