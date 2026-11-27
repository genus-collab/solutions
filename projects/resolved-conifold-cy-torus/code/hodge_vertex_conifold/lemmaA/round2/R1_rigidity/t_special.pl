#!/usr/bin/perl
# t_special.pl -- the special-framing lemma: at tau0 = -k0/q (gcd(k0,q) = 1, 1 <= k0 < q), genus G, the first-order two-leaf
# regularity conditions  Y(m,n)(tau0) = 0 for q | n, m >= 1, m + n <= 2G-1  (Y(m,n) = [u1^m u2^n] y1 y2 (F(y1)-F(y2))/(y1-y2),
# F(y) = sum_c c y_{c-1} u(y)^c, u = y(1-y)^tau0) have the same solution space in (y_1..y_{2G-3}) as
#   [y^N] (1-y)^(k0-1) (1-(1+tau0) y) F(y) = 0  for q <= N <= 2G-2.
# Exact (Math::BigRat).  Usage: perl t_special.pl G   (all q <= 2G-2, k0 < q coprime)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($G) = @ARGV; my $L = 2*$G - 3; my $Kmax = 2*$G - 2;
sub gbin { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
sub rank { my @M = map { [ @$_ ] } @{$_[0]}; my $nc = $_[1]; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { unless ($M[$i][$c]->is_zero) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = R(1)/$M[$rk][$c]; $M[$rk] = [ map { $_ * $inv } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || $M[$i][$c]->is_zero; my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = $M[$i][$j] - $f*$M[$rk][$j]; } } $rk++; } return $rk; }
my ($tests, $bad) = (0, 0);
for my $q (2..$Kmax) { for my $k0 (1..$q-1) { next unless gcd($k0, $q) == 1; my $t0 = R(-$k0)/R($q);
  # f_k as row vectors over y_1..y_L (index c-1 -> column c-2)
  my @f; for my $k (1..$Kmax) { my @r = map { R(0) } 1..$L; for my $c (2..$k) { next if $c - 1 > $L; $r[$c-2] = R($c) * R((-1)**($k-$c)) * gbin($c*$t0, $k-$c); } $f[$k] = \@r; }
  my %cc; for my $n (1..$Kmax) { for my $i (1..$n) { $cc{"$n,$i"} = R($i)/R($n) * gbin($t0*$n + $n - $i - 1, $n - $i); } }
  my @rows1; for (my $n = $q; $n <= $Kmax; $n += $q) { for my $m (1..$Kmax-$n+1) { next if $m + $n > 2*$G - 1; my @r = map { R(0) } 1..$L;
      for my $j (1..$m) { for my $i (1..$n) { my $k = $i + $j - 1; my $w = $cc{"$m,$j"} * $cc{"$n,$i"}; next if $w->is_zero; $r[$_] = $r[$_] + $w*$f[$k][$_] for 0..$L-1; } }
      push @rows1, \@r; } }
  # p-conditions: (1-y)^(k0-1)(1-(1+t0)y) as a polynomial
  my @w = (R(1)); for (1..$k0-1+($ENV{EXP_SHIFT}//0)) { my @nw = map { R(0) } 0..@w; for my $i (0..$#w) { $nw[$i] = $nw[$i] + $w[$i]; $nw[$i+1] = $nw[$i+1] - $w[$i]; } @w = @nw; }
  { my @nw = map { R(0) } 0..@w; for my $i (0..$#w) { $nw[$i] = $nw[$i] + $w[$i]; $nw[$i+1] = $nw[$i+1] - (1+$t0)*$w[$i]; } @w = @nw; }
  my @rows2; for my $N ($q..$Kmax) { my @r = map { R(0) } 1..$L; for my $i (0..$#w) { my $k = $N - $i; next if $k < 1; $r[$_] = $r[$_] + $w[$i]*$f[$k][$_] for 0..$L-1; } push @rows2, \@r; }
  my $r1 = @rows1 ? rank(\@rows1, $L) : 0; my $r2 = @rows2 ? rank(\@rows2, $L) : 0; my $r12 = rank([ @rows1, @rows2 ], $L);
  $tests++; my $ok = ($r1 == $r2 && $r2 == $r12); $bad++ unless $ok;
  printf "G=%d q=%d k0=%d: rank two-leaf %d, rank p-conditions %d, joint %d %s\n", $G, $q, $k0, $r1, $r2, $r12, ($ok ? "SAME" : "DIFFERENT") if !$ok || $ENV{VERBOSE}; } }
print "t_special G=$G: $tests framings, $bad differ\n";
