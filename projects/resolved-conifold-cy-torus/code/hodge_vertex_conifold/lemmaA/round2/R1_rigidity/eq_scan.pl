#!/usr/bin/perl
# eq_scan.pl -- for every a <= amax, k0 coprime to a (k0 <= a/2), odd D with D+2 <= a: smallest G at which the closed-form forcing
# matrix M_{N,n} (eq_kernel.pl) forces p_{D+1}; 'never' if not forced by the exact functional equation (N up to 4a+20).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($amax, $Dmax) = @ARGV; $amax //= 20; $Dmax //= 9;
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
sub gb { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
sub rank { my ($M0, $nc) = @_; my @M = map { [ @$_ ] } @$M0; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { unless ($M[$i][$c]->is_zero) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = R(1)/$M[$rk][$c]; $M[$rk] = [ map { $_ * $inv } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || $M[$i][$c]->is_zero; my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = $M[$i][$j] - $f*$M[$rk][$j]; } } $rk++; } return $rk; }
for (my $D = 3; $D <= $Dmax; $D += 2) { my @line;
  for my $a ($D+2..$amax) { for my $k0 (1..int($a/2)) { next unless gcd($k0, $a) == 1;
      my @cols = ($D+1 .. $a-1); my $nc = @cols; my @rows; my $found = 'never';
      for (my $N = $D + 2; $N <= 4*$a + 20; $N += 2) { my @r = map { my $n = $_; $N < $n ? R(0) : R($a-$k0) * gb($k0 + $N - $n - 1, $N - $n) + R($k0) * gb($k0 - 1 - $n, $N - $n) } @cols; push @rows, \@r;
        my $r1 = rank(\@rows, $nc); my @rows2 = map { [ @{$_}[1..$nc-1] ] } @rows; my $r2 = $nc > 1 ? rank(\@rows2, $nc-1) : 0;
        if ($r1 - $r2 == 1) { $found = ($N + 3)/2; last; } }
      my $g0 = $a - ($D - 1)/2; push @line, sprintf("%d/%d:%s%s", $k0, $a, $found, ($found ne 'never' && $found != $g0) ? "(*)" : ""); } }
  print "D=$D (threshold a-(D-1)/2 marked (*) if different): @line\n"; }
