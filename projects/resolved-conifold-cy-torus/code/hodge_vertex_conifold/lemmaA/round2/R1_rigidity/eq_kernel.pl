#!/usr/bin/perl
# eq_kernel.pl -- the pointwise forcing problem at tau0 = -k0/a in closed form.
# Two-leaf regularity at tau0:  F(y; tau0) = p(y)/((1-y)^{k0-1}(1-(1+tau0)y)) mod y^{2G-1}, deg p <= a-1, p = O(y^{D+1}).
# Closed-form relations (closed form of the polar-order relations):  the odd part of
#     a R(s) = (a-k0)(1-s)^{-k0} p(s) + k0 (1+s)^{k0-1} p(s/(1+s))
# vanishes mod s^{2G-1}, i.e. for every odd N <= 2G-3:  sum_n M_{N,n} p_n = 0,
#     M_{N,n} = (a-k0) binom(k0+N-n-1, N-n) + k0 binom(k0-1-n, N-n)    (generalised binomial).
# p_{D+1} = F_{D+1}(tau0) = (D+1) y_D(tau0).  Forcing <=> column n = D+1 not in the span of the other columns.
# Usage: perl eq_kernel.pl a k0 D G   (G = 0: exact functional equation, N up to a large bound)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($a, $k0, $D, $G) = @ARGV;
my $Nmax = $G ? 2*$G - 3 : 6*$a + 20;
sub gb { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
my @cols = ($D+1 .. $a-1);
my @rows;
for (my $N = 1; $N <= $Nmax; $N += 2) { my @r = map { my $n = $_; $N < $n ? R(0) : R($a-$k0) * gb($k0 + $N - $n - 1, $N - $n) + R($k0) * gb($k0 - 1 - $n, $N - $n) } @cols; push @rows, \@r if grep { !$_->is_zero } @r; }
sub rank { my @M = map { [ @$_ ] } @{$_[0]}; my $nc = $_[1]; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { unless ($M[$i][$c]->is_zero) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = R(1)/$M[$rk][$c]; $M[$rk] = [ map { $_ * $inv } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || $M[$i][$c]->is_zero; my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = $M[$i][$j] - $f*$M[$rk][$j]; } }
    $rk++; } return $rk; }
my $nc = @cols; my $r1 = rank(\@rows, $nc); my @rows2 = map { [ @{$_}[1..$nc-1] ] } @rows; my $r2 = $nc > 1 ? rank(\@rows2, $nc - 1) : 0;
printf "a=%d k0=%d D=%d G=%s: unknowns p_%d..p_%d (%d), odd-order equations %d, rank %d, kernel dim %d, p_{D+1} forced: %s\n",
  $a, $k0, $D, ($G || 'inf'), $D+1, $a-1, $nc, scalar(@rows), $r1, $nc - $r1, ($r1 - $r2 == 1 ? 'YES' : 'no');
