#!/usr/bin/perl
# t_T1_modp.pl -- modular test of the closed form (T1) of the polar-order relations to high order:
#   relation_k = sum_L beta_{k,L} y_L  with  beta_{k,L} = kappa_k (-1)^d [ (1+tau) binom(N tau, d) - tau binom(N(1+tau)+d-1, d) ],
#   N = L+1, d = k-L, kappa_1 = tau + 1/2, kappa_k = B_k/(k k!) (k >= 2; zero for odd k >= 3).
# The left side is computed by kpmod.pm (full KP propagation + leaf curve), mod p, at the given framing values.
# Usage: perl t_T1_modp.pl K tau1 tau2 ...   (taus as fractions a/b, may be special framings)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod;
my ($K, @taus) = @ARGV;
# Bernoulli numbers mod p
my @fct = (1); $fct[$_] = mm($fct[$_-1], $_) for 1..$K+2;
sub binm { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; return mm($fct[$n], mi(mm($fct[$k], $fct[$n-$k]))); }
my @B = (1); for my $m (1..$K) { my $s = 0; for my $k (0..$m-1) { $s = ($s + mm(binm($m+1, $k), $B[$k])) % $P; } $B[$m] = md(-mm($s, mi($m+1))); }
my ($n, $bad) = (0, 0);
for my $ts (@taus) { my ($a, $b) = split m{/}, $ts; $b //= 1; my $tau = mfrac($a, $b);
  my $S = kp_setup($tau, 1, $K); my $rows = kp_rows($S);
  for my $r (@$rows) { my ($k, $eq) = @$r;
    my $kap = $k == 1 ? md($tau + mi(2)) : mm($B[$k], mi(mm($k, $fct[$k])));
    for my $L (1..$k) { my $N = $L + 1; my $d = $k - $L;
      # generalised binomials binom(x, d) = x(x-1)...(x-d+1)/d!
      my $g1 = 1; my $x1 = mm($N, $tau); for my $i (0..$d-1) { $g1 = mm($g1, md($x1 - $i)); } $g1 = mm($g1, mi($fct[$d]));
      my $g2 = 1; my $x2 = md(mm($N, md(1 + $tau)) + $d - 1); for my $i (0..$d-1) { $g2 = mm($g2, md($x2 - $i)); } $g2 = mm($g2, mi($fct[$d]));
      my $pred = mm($kap, mm(md(($d % 2) ? -1 : 1), md(mm(md(1 + $tau), $g1) - mm($tau, $g2))));
      $n++; if ($pred != $eq->[$L-1]) { $bad++; print "MISMATCH tau=$ts k=$k L=$L\n" if $bad < 20; } } }
  print "tau=$ts: checked orders 1..$K\n"; }
print "t_T1_modp: $n coefficients compared, $bad mismatches\n";
