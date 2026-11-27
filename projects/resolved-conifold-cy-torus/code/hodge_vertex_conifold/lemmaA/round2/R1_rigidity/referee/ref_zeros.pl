#!/usr/bin/perl
# ref_zeros.pl -- referee's recount of the zero count (route s.10).
# (1) Z1(D): zeros (with multiplicity, in sigma) of sigma^G Xi_G(D+1) forced by Corollary E (merge formula), computed by
#     brute force over all T at level D for D <= 14: at tau0 = -c/a the zero order in tau is
#     max_T [ #{t in T : a | t} ] - [a | D+1]; at a = 2 (sigma0 = -1/4) the sigma-order is ceil(tau-order/2);
#     compared with the closed formula Z1(D) = floor(D/2) + Sum_{a=3}^{D+1} (phi(a)/2) max(0, floor(D/(a-1)) - [a|D+1]).
# (2) Z2(G,D) = #{(k0,q) : gcd = 1, 1 <= k0 <= D+1, 2k0 < q, q >= D+2, 3(q-k0) < 2G}: direct enumeration (G <= 60)
#     and divisor (Moebius) counting (G <= Gmax).  Report every (G, odd D in [3, 2G-3]) with Z1 + Z2 < G-1.
# (3) The analytic bounds (B1), (B2) checked numerically.
use strict; use warnings;
my ($Gmax) = @ARGV; $Gmax //= 400;
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
my $NMAX = 2*$Gmax + 10;
my @phi = (0, 1); for my $n (2..$NMAX) { my $c = 0; my $m = $n; my $r = $n; for (my $p = 2; $p*$p <= $m; $p++) { if ($m % $p == 0) { $r = $r / $p * ($p - 1); $m /= $p while $m % $p == 0; } } $r = $r / $m * ($m - 1) if $m > 1; $phi[$n] = $r; }
sub Z1f { my $D = shift; my $z = int($D/2); for my $a (3..$D+1) { my $m = int($D/($a-1)) - ((($D+1) % $a) == 0 ? 1 : 0); $z += ($phi[$a]/2) * $m if $m > 0; } return $z; }
# brute force over partitions
sub parts { my ($n, $max) = @_; return ([]) if $n == 0; my @o; for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for parts($n - $k, $k); } return @o; }
my $bfbad = 0;
for my $D (1..14) { my @Ts = map { [ map { $_ + 1 } @$_ ] } parts($D, $D);
  my $z = 0;
  for my $a (2..$D+1) { my $best = 0; for my $T (@Ts) { my $c = grep { $_ % $a == 0 } @$T; $best = $c if $c > $best; }
    my $ord = $best - ((($D+1) % $a) == 0 ? 1 : 0); next if $ord <= 0;
    if ($a == 2) { $z += int(($ord + 1)/2); } else { $z += ($phi[$a]/2) * $ord; } }
  my $f = Z1f($D); $bfbad++ if $f != $z; printf "D=%2d Z1 brute %3d formula %3d\n", $D, $z, $f if $D <= 10 || $f != $z; }
print "Z1 brute force vs formula (D <= 14): $bfbad differences\n";
# divisors
my @divs; for my $n (1..$NMAX) { for (my $m = $n; $m <= $NMAX; $m += $n) { push @{$divs[$m]}, $n; } }
my @mu = (0, 1); for my $n (2..$NMAX) { my $m = $n; my $k = 0; my $sq = 0; for (my $p = 2; $p*$p <= $m; $p++) { if ($m % $p == 0) { $m /= $p; $k++; if ($m % $p == 0) { $sq = 1; last; } } } $k++ if !$sq && $m > 1; $mu[$n] = $sq ? 0 : ($k % 2 ? -1 : 1); }
# recompute mu properly (the loop above is sloppy for repeated small primes); use a sieve instead
@mu = (0) x ($NMAX+1); $mu[1] = 1; for my $n (1..$NMAX) { next unless $mu[$n]; for (my $m = 2*$n; $m <= $NMAX; $m += $n) { $mu[$m] -= $mu[$n]; } }
sub coprime_count { my ($k, $lo, $hi) = @_; return 0 if $hi < $lo; my $s = 0; for my $d (@{$divs[$k]}) { next unless $mu[$d]; $s += $mu[$d] * (int($hi/$d) - int(($lo-1)/$d)); } return $s; }
sub Z2m { my ($G, $D) = @_; my $z = 0; for my $k0 (1..$D+1) { my $lo = $D + 2; $lo = 2*$k0 + 1 if 2*$k0 + 1 > $lo; my $hi = $k0 + int((2*$G - 1)/3);   # 3(q-k0) < 2G <=> q - k0 <= floor((2G-1)/3)
    $z += coprime_count($k0, $lo, $hi); } return $z; }
sub Z2d { my ($G, $D) = @_; my $z = 0; for my $q ($D+2 .. 2*$G) { for my $k0 (1..$D+1) { next unless 2*$k0 < $q && 3*($q - $k0) < 2*$G && gcd($k0, $q) == 1; $z++; } } return $z; }
my ($fails, $minmarg, $argmin) = (0, 1e9, '');
my $dbad = 0;
for my $G (3..$Gmax) { for (my $D = 3; $D <= 2*$G - 3; $D += 2) {
    my $z2 = Z2m($G, $D); if ($G <= 60) { my $zd = Z2d($G, $D); $dbad++ if $zd != $z2; }
    my $tot = Z1f($D) + $z2; my $marg = $tot - ($G - 1);
    if ($marg < 0) { $fails++; print "count fails: G=$G D=$D Z1=", Z1f($D), " Z2=$z2 (need ", $G-1, ")\n"; }
    elsif ($G >= 8 && $marg < $minmarg) { $minmarg = $marg; $argmin = "G=$G D=$D"; }
    # (B1), (B2)
    my $b1 = (16/9)*$G - (8/3)*($D+1) - 5/3; if ($z2 < $b1 - 1e-9) { print "B1 violated G=$G D=$D: Z2=$z2 < $b1\n"; }
    if ($D >= 8) { my $b2 = 0.0376*$D*$D + 0.5*$D - 2.75; if (Z1f($D) < $b2 - 1e-9) { print "B2 violated D=$D\n"; } } } }
print "Z2 Moebius vs direct (G <= 60): $dbad differences\n";
print "zero count G=3..$Gmax: $fails failures (all listed above); minimal margin for G >= 8: $minmarg at $argmin\n";
