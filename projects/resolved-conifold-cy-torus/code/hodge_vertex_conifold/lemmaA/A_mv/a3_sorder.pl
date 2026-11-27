#!/usr/bin/perl
# a3_sorder.pl -- structure of Q_h(s;R) = Sum_k s^k/k! <Omega lam_R Tt^k> (c = 1) as a polynomial in s:
# order of vanishing at s = -1, and the values of the derivatives at s = -1 (Taylor coefficients in (1+s)).
# Usage: perl a3_sorder.pl h "R1;R2;..." [taus]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($h, $Rlist, $taulist) = @ARGV;
my @Rs = map { [ split /,/, $_ ] } split /;/, $Rlist;
my @taus = map { R($_) } split /,/, ($taulist // '2/7,3/2');
for my $tau (@taus) {
  for my $Rp (@Rs) { my $r = scalar @$Rp; my $N = 3*$h - 3 + $r; my $J = 0; $J += $_ for @$Rp; my $D = $J - $r;
    my $M = mv_unit($Rp, $tau, $h, $N);
    my @q; for my $m (0..$N) { my $s = R(0); for my $a (0..$m) { $s += $M->[$h][$a] / fact($a) * R((-1)**($m-$a)) * lag($m-$a, 2*$h-3+$r+$a, -($J+$a)); } $q[$m] = $s; }
    # Taylor coefficients at s = -1: Q(s) = sum_k t_k (1+s)^k,  t_k = sum_m q_m binom(m,k) (-1)^{m-k}
    my @t; for my $k (0..$N) { my $s = R(0); for my $m ($k..$N) { $s += $q[$m] * binom($m, $k) * R((-1)**($m-$k)); } $t[$k] = $s; }
    my $ord = 0; $ord++ while $ord <= $N && $t[$ord]->is_zero;
    printf "h=%d tau=%s R=(%s) D=%d  2h-2-D=%d : ord_{s=-1} Q = %s ;  first Taylor coeffs: %s\n", $h, $tau, join(',', @$Rp), $D, 2*$h-2-$D,
      ($ord > $N ? 'inf' : $ord), join(', ', map { $t[$_] } 0..($N < 5 ? $N : 5)); } }
