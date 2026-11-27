#!/usr/bin/perl
# count_zeros.pl -- the zero count that closes (OL) at level D (odd, 3 <= D <= 2G-3), genus G, given the forcing theorem.
#   Z1(D) = Route A's merge-formula count (zeros with multiplicity of sigma^G Xi_G(D+1) at the framings -c/a with a <= D+1):
#           Z1 = floor(D/2) + sum_{a=3}^{D+1} (phi(a)/2) * max(0, floor(D/(a-1)) - [a | D+1])
#   Z2(G,D) = #{ reduced k0/a : 1 <= k0 <= D+1, 2 k0 < a, a >= D+2, 3(a - k0) < 2G }   (simple zeros, forcing theorem)
# (OL) at (G, D) follows if Z1 + Z2 >= G - 1 (sigma^G Xi_G(D+1) = sigma Q, deg Q <= G-2).
# Usage: perl count_zeros.pl Gmin Gmax      prints the minimal margin Z1 + Z2 - (G-1) over D for each G, and any failure.
use strict; use warnings;
my ($Gmin, $Gmax) = @ARGV; $Gmin //= 3; $Gmax //= 200;
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
my @phi = (0, 1); for my $n (2..4*$Gmax) { my $c = 0; my $m = $n; my $r = $n; for (my $p = 2; $p*$p <= $m; $p++) { if ($m % $p == 0) { $r = $r / $p * ($p - 1); $m /= $p while $m % $p == 0; } } $r = $r / $m * ($m - 1) if $m > 1; $phi[$n] = $r; }
my %Z1c;
sub Z1 { my $D = shift; return $Z1c{$D} //= do { my $z = int($D/2); for my $a (3..$D+1) { my $o = int($D/($a-1)) - ((($D+1) % $a) == 0 ? 1 : 0); $z += ($phi[$a]/2) * $o if $o > 0; } $z }; }
my $fail = 0;
for my $G ($Gmin..$Gmax) { my ($minm, $argD) = (1e18, 0);
  for (my $D = 3; $D <= 2*$G - 3; $D += 2) {
    my $z2 = 0; my $amax = int((2*$G - 1)/3) + $D + 1;   # 3(a-k0) < 2G with k0 <= D+1
    for my $a ($D+2 .. $amax) { for my $k0 (1 .. ($D+1 < ($a-1)/2 ? $D+1 : int(($a-1)/2))) { next unless 3*($a - $k0) < 2*$G; next unless gcd($k0, $a) == 1; $z2++; } }
    my $m = Z1($D) + $z2 - ($G - 1); if ($m < $minm) { ($minm, $argD) = ($m, $D); }
    if ($m < 0) { $fail++; print "FAIL G=$G D=$D: Z1=", Z1($D), " Z2=$z2 need ", $G-1, "\n"; } }
  print "G=$G: minimal margin $minm (at D=$argD)\n" if $G <= 40 || $G % 50 == 0; }
print "count_zeros: G=$Gmin..$Gmax, failures: $fail\n";
