#!/usr/bin/perl
# beta_to_f.pl -- closed-form relation rows in the basis f_m = [X^m] V, X = 1/lambda = y/(1 + tau y)  (V(z) = sum_c c y_{c-1} u^c, u = 1/z).
# c y_{c-1} = [u^c] V = sum_m f_m [u^c] X(u)^m; computed via y: X = y/(1+tau y), u = y(1-y)^tau.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($file, $taus) = @ARGV; my $tau = R($taus);
my @rows; open my $fh, '<', $file or die; while (<$fh>) { if (/^REL (\d+): (.*)$/) { $rows[$1] = [ map { R($_) } split ' ', $2 ]; } } close $fh;
my $K = $#rows; my $N = $K + 2;
# series of y(u): inverse of u = y(1-y)^tau
my $uy = series_mul([R(0), R(1)], series_pow([R(1), R(-1)], $tau, $N), $N);   # u(y)
# invert: y(u) by iteration y = u (1-y)^{-tau}
my @y = (R(0), R(1), map { R(0) } 2..$N);
for (1..$N) { my @omy = map { -$_ } @y; $omy[0] = R(1); my $p = series_pow(\@omy, -$tau, $N); my @ny = (R(0), @$p[0..$N-1]); @y = @ny; }
my @opty = map { $_ * $tau } @y; $opty[0] = R(1);
my $X = series_mul(\@y, series_inv(\@opty, $N), $N);   # X(u)
my @Xp = ([R(1), map { R(0) } 1..$N]); for my $m (1..$N) { $Xp[$m] = series_mul($Xp[$m-1], $X, $N); }
for my $k (1..$K) { my @g;
  for my $m (1..$K+1) { my $s = R(0);
    for my $c ($m..$K+1) { my $L = $c - 1; next if $L < 1; my $b = $rows[$k][$L-1]; next if $b->is_zero; $s = $s + $b * $Xp[$m][$c] / $c; }
    push @g, $s; }
  print "k=$k: ", join("  ", map { "f$_:" . $g[$_-1] } grep { !$g[$_-1]->is_zero } 1..$K+1), "\n"; }
