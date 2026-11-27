#!/usr/bin/perl
# ref_zeros.pl -- referee 2: independent zero count for s.10.
# (1) Z1(D) recomputed from first principles: for every reduced c/a (2 <= a <= D+1) the forced vanishing order of Xi_G(D+1)
#     at tau = -c/a is max over ALL multisets T (parts >= 2, level D) of ord(prod Pi_t) - ord(Pi_{D+1}) (orders of the linear
#     factors j tau + b at tau = -c/a), converted to sigma-order (halved at a = 2), one sigma-point per pair {c, a-c}.
#     Compared with the author's closed formula floor(D/2) + sum_a (phi(a)/2) max(0, floor(D/(a-1)) - [a | D+1]).
# (2) Z2(G,D) by direct enumeration of reduced k0/q with k0 <= D+1, 2k0 < q, q >= D+2, 3(q-k0) < 2G.
# (3) margin Z1 + Z2 - (G-1) for odd D in [3, 2G-3], G in [Gmin, Gmax]; reports failures and the minimal margin.
# Usage: perl ref_zeros.pl Gmin Gmax [Dbrute]   (Dbrute: largest D for the brute-force Z1 check, default 14)
use strict; use warnings;
my ($Gmin, $Gmax, $Db) = @ARGV; $Db //= 14;
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
my $NM = 2*$Gmax + 4; my @phi = (0, 1); for my $n (2..$NM) { my $c = 0; for my $k (1..$n) { $c++ if gcd($k, $n) == 1; } $phi[$n] = $c; }
sub Z1f { my $D = shift; my $z = int($D/2); for my $a (3..$D+1) { my $o = int($D/($a-1)) - (($D+1) % $a == 0 ? 1 : 0); $z += ($phi[$a]/2)*$o if $o > 0; } return $z; }
# brute force Z1
sub parts { my ($n, $mx) = @_; return ([]) if $n == 0; my @o; for my $k (reverse 1..($mx < $n ? $mx : $n)) { push @o, [$k, @$_] for parts($n - $k, $k); } return @o; }
sub ordPi { my ($j, $c, $a) = @_; my $o = 0; for my $b (1..$j-1) { $o++ if $j*$c == $b*$a; } return $o; }   # zeros of j tau + b at tau = -c/a
my $bad1 = 0;
for my $D (1..$Db) { my @Ts = map { [ map { $_ + 1 } @$_ ] } parts($D, $D); my $z = 0;
  for my $a (2..$D+1) { for my $c (1..int($a/2)) { next unless gcd($c, $a) == 1; next if 2*$c == $a && $a != 2;
      my $best = 0; for my $T (@Ts) { my $o = 0; $o += ordPi($_, $c, $a) for @$T; $o -= ordPi($D+1, $c, $a); $best = $o if $o > $best; }
      $z += ($a == 2) ? int($best/2) : $best; } }
  my $zf = Z1f($D); if ($z != $zf) { $bad1++; print "Z1 MISMATCH D=$D: brute $z, formula $zf\n"; } else { print "D=$D: Z1 = $z\n" if $D <= 12; } }
print "brute-force Z1 check D <= $Db: $bad1 mismatches\n";
my $fail = 0; my %minm;
for my $G ($Gmin..$Gmax) { my ($mm, $mD) = (1e9, 0);
  for (my $D = 3; $D <= 2*$G-3; $D += 2) { my $z2 = 0;
    for my $k0 (1..$D+1) { for my $q (($D+2 > 2*$k0+1 ? $D+2 : 2*$k0+1)..($k0 + int((2*$G - 1)/3))) { next unless 3*($q - $k0) < 2*$G; $z2++ if gcd($k0, $q) == 1; } }
    my $m = Z1f($D) + $z2 - ($G - 1); ($mm, $mD) = ($m, $D) if $m < $mm;
    if ($m < 0) { $fail++; print "FAIL G=$G D=$D: Z1=", Z1f($D), " Z2=$z2 need ", $G-1, "\n"; } }
  print "G=$G: min margin $mm at D=$mD\n" if $G <= 24 || $G % 50 == 0; }
print "ref_zeros: G=$Gmin..$Gmax, $fail failures\n";
