#!/usr/bin/perl
# count_zeros2.pl -- fast version of count_zeros.pl (same counts):
#   Z1(D) = floor(D/2) + sum_{a=3}^{D+1} (phi(a)/2) max(0, floor(D/(a-1)) - [a | D+1])            (Route A, Cor. E)
#   Z2(G,D) = sum_{k0=1}^{D+1} #{ q : gcd(q,k0) = 1, max(D+2, 2k0+1) <= q <= k0 + ceil(2G/3) - 1 }   (forcing theorem)
# (3(q-k0) < 2G  <=>  q - k0 <= ceil(2G/3) - 1.)  Counts coprime integers by Moebius inversion over the divisors of k0.
# Prints min margin Z1 + Z2 - (G-1) over odd D in [3, 2G-3] per G; also checks the analytic bounds used for large G:
#   (B1) Z2(G,D) >= (16/9) G - (8/3)(D+1) - 5/3 whenever ceil(2G/3) - 1 >= D + 1,
#   (B2) Z1(D) >= 0.037595 D^2 + 0.5 D - 2.75 for D >= 8.
use strict; use warnings;
my ($Gmin, $Gmax) = @ARGV; $Gmin //= 3; $Gmax //= 300;
my $NM = 4*$Gmax + 10;
my @phi = (0) x ($NM+1); my @mu = (0) x ($NM+1); my @divs; 
{ my @lp = (0) x ($NM+1); my @pr; $mu[1] = 1; $phi[1] = 1;
  for my $i (2..$NM) { if (!$lp[$i]) { $lp[$i] = $i; push @pr, $i; $phi[$i] = $i - 1; $mu[$i] = -1; }
    for my $p (@pr) { last if $p > $lp[$i] || $p*$i > $NM; $lp[$p*$i] = $p; if ($p == $lp[$i]) { $phi[$p*$i] = $phi[$i]*$p; $mu[$p*$i] = 0; } else { $phi[$p*$i] = $phi[$i]*($p-1); $mu[$p*$i] = -$mu[$i]; } } } }
for my $d (1..$NM) { next unless $mu[$d]; for (my $m = $d; $m <= $NM; $m += $d) { push @{$divs[$m]}, $d; } }
sub cop { my ($k, $x) = @_; return 0 if $x <= 0; my $s = 0; $s += $mu[$_] * int($x/$_) for @{$divs[$k]}; return $s; }   # #{1<=q<=x : gcd(q,k)=1}
my @Z1; for my $D (1..2*$Gmax) { my $z = int($D/2); for my $a (3..$D+1) { my $o = int($D/($a-1)) - ((($D+1) % $a) == 0 ? 1 : 0); $z += ($phi[$a]/2) * $o if $o > 0; } $Z1[$D] = $z; }
my ($fail, $b1fail, $b2fail) = (0, 0, 0);
for my $D (8..2*$Gmax) { $b2fail++ if $Z1[$D] < 0.037595*$D*$D + 0.5*$D - 2.75; }
for my $G ($Gmin..$Gmax) { my ($minm, $argD) = (1e18, 0); my $g = int((2*$G + 2)/3) - 1;   # ceil(2G/3) - 1
  for (my $D = 3; $D <= 2*$G - 3; $D += 2) {
    my $z2 = 0; for my $k0 (1..$D+1) { my $lo = $D + 2; $lo = 2*$k0 + 1 if 2*$k0 + 1 > $lo; my $hi = $k0 + $g; next if $hi < $lo; $z2 += cop($k0, $hi) - cop($k0, $lo - 1); }
    $b1fail++ if $g >= $D + 1 && $z2 < (16/9)*$G - (8/3)*($D+1) - 5/3;
    my $m = $Z1[$D] + $z2 - ($G - 1); if ($m < $minm) { ($minm, $argD) = ($m, $D); }
    if ($m < 0) { $fail++; print "FAIL G=$G D=$D: Z1=$Z1[$D] Z2=$z2 need ", $G-1, "\n"; } }
  print "G=$G: minimal margin $minm (at D=$argD)\n" if $G <= 30 || $G % 100 == 0; }
print "count_zeros2: G=$Gmin..$Gmax, failures: $fail; bound (B1) violations: $b1fail; bound (B2) violations (D=8..", 2*$Gmax, "): $b2fail\n";
