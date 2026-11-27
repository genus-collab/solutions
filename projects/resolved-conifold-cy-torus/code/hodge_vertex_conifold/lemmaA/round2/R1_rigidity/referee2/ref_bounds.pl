#!/usr/bin/perl
# ref_bounds.pl -- referee 2: numerical check of the analytic bounds (B1), (B2) of s.10 and of the case split for large G.
#   (B1) Z2(G,D) >= (16/9)G - (8/3)(D+1) - 5/3 ;  (B2) Z1(D) >= 0.0376 D^2 + 0.5 D - 2.75 (D >= 8)
use strict; use warnings;
my ($Dmax, $Gmax) = @ARGV; my $NM = $Dmax + 2;
my @phi = (0..$NM); for my $p (2..$NM) { next unless $phi[$p] == $p; for (my $m = $p; $m <= $NM; $m += $p) { $phi[$m] -= $phi[$m]/$p; } }
my $b2 = 0; my @Z1;
for my $D (1..$Dmax) { my $z = int($D/2); for my $a (3..$D+1) { my $o = int($D/($a-1)) - (($D+1) % $a == 0 ? 1 : 0); $z += ($phi[$a]/2)*$o if $o > 0; } $Z1[$D] = $z; $b2++ if $D >= 8 && $z < 0.0376*$D*$D + 0.5*$D - 2.75; }
print "(B2) for 8 <= D <= $Dmax: $b2 violations\n";
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
my ($b1, $fail) = (0, 0);
for my $G (map { int($_) } (300, 301, 500, 777, 1001, 1500, $Gmax)) { my $g = int((2*$G + 2)/3) - 1;
  for (my $D = 3; $D <= 2*$G - 3; $D += 2) { last if $D > $Dmax;
    my $z2 = 0; for my $k0 (1..$D+1) { my $lo = $D + 2; $lo = 2*$k0 + 1 if 2*$k0 + 1 > $lo; for my $q ($lo..$k0 + $g) { $z2++ if gcd($k0, $q) == 1; } }
    $b1++ if $z2 < (16/9)*$G - (8/3)*($D+1) - 5/3; $fail++ if $Z1[$D] + $z2 < $G - 1; } }
print "(B1) and total count at sample G (300..$Gmax, all odd D): $b1 (B1) violations, $fail count failures\n";
