#!/usr/bin/perl
# ref_margin.pl -- referee 2: minimal margin m(G) = min over odd D in [3, 2G-3] of Z1(D) + Z2(G,D) - (G-1), for G <= Gmax,
# with Z2 counted by inclusion-exclusion over the prime divisors of k0 (own implementation). Reports G with m(G) < 0 and m(G) < 1.
use strict; use warnings;
my ($Gmax) = @ARGV; my $NM = 2*$Gmax + 4;
my @phi = (0..$NM); my @pf = map { [] } 0..$NM;
for my $p (2..$NM) { next unless $phi[$p] == $p; for (my $m = $p; $m <= $NM; $m += $p) { $phi[$m] -= $phi[$m]/$p; push @{$pf[$m]}, $p; } }
sub cop { my ($k, $x) = @_; return 0 if $x <= 0; my @p = @{$pf[$k]}; my $s = 0; for my $mask (0..(1 << @p) - 1) { my ($d, $b) = (1, 0); for my $i (0..$#p) { if ($mask >> $i & 1) { $d *= $p[$i]; $b++; } } $s += (($b % 2) ? -1 : 1) * int($x/$d); } return $s; }
my @Z1; for my $D (1..2*$Gmax) { my $z = int($D/2); for my $a (3..$D+1) { my $o = int($D/($a-1)) - (($D+1) % $a == 0 ? 1 : 0); $z += ($phi[$a]/2)*$o if $o > 0; } $Z1[$D] = $z; }
my (@neg, @zero); my $minG8 = 1e9;
for my $G (3..$Gmax) { my $g = int((2*$G + 2)/3) - 1; my $mm = 1e9;
  for (my $D = 3; $D <= 2*$G - 3; $D += 2) { my $z2 = 0; for my $k0 (1..$D+1) { my $lo = $D + 2; $lo = 2*$k0 + 1 if 2*$k0 + 1 > $lo; my $hi = $k0 + $g; next if $hi < $lo; $z2 += cop($k0, $hi) - cop($k0, $lo - 1); }
    my $m = $Z1[$D] + $z2 - ($G - 1); $mm = $m if $m < $mm; }
  push @neg, $G if $mm < 0; push @zero, $G if $mm == 0; $minG8 = $mm if $G >= 8 && $mm < $minG8; }
print "ref_margin G<=$Gmax: m(G) < 0 for G in {@neg}; m(G) = 0 for G in {@zero}; min over 8 <= G <= $Gmax of m(G): $minG8\n";
