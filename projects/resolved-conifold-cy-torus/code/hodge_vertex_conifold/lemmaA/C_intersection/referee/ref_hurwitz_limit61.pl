#!/usr/bin/perl
# Byproduct check (not part of Theorem H): the conjectural Hurwitz limit of (6.1), proof_M_residue.md s.6a(iii),
#   Xi^H_h(j) = prod_{l=1}^{2h-2} (j-l) * ((j-2h)/j) * [x^{2h}] H(x)^j,  log H = Sum_k (2k-1) B_{2k} x^{2k} / (2k (2k)!),
# evaluated for comparison with the TR-free Hurwitz-number values of ref_elsv.pl.  Usage: perl ref_hurwitz_limit61.pl h j
use strict; use warnings; use Math::BigRat; use Math::BigInt;
my ($h, $j) = @ARGV; my $n = 2 * $h;
# Bernoulli numbers B_0..B_n
my @B = (Math::BigRat->new(1)); for my $m (1 .. $n) { my $s = Math::BigRat->new(0); for my $k (0 .. $m - 1) { my $bin = Math::BigInt->new($m + 1)->bnok($k); $s += Math::BigRat->new($bin->bstr) * $B[$k]; } $B[$m] = -$s / ($m + 1); }
sub fct { Math::BigInt->new(shift)->bfac() }
my @L = map { Math::BigRat->new(0) } 0 .. $n; for my $k (1 .. $h) { $L[2*$k] = Math::BigRat->new(2*$k - 1) * $B[2*$k] / (Math::BigRat->new(2*$k) * Math::BigRat->new(fct(2*$k)->bstr)); }
# H^j = exp(j L)
my @E = (Math::BigRat->new(1), map { Math::BigRat->new(0) } 1 .. $n); my @term = @E;
for my $p (1 .. $n) { my @t = map { Math::BigRat->new(0) } 0 .. $n; for my $a (0 .. $n) { next if $term[$a]->is_zero; for my $b (0 .. $n - $a) { next if $L[$b]->is_zero; $t[$a + $b] += $term[$a] * $L[$b] * $j / $p; } } @term = @t; $E[$_] += $term[$_] for 0 .. $n; }
my $pre = Math::BigRat->new(1); $pre *= ($j - $_) for 1 .. 2 * $h - 2; $pre *= Math::BigRat->new($j - 2 * $h) / $j;
print "(6.1)-Hurwitz limit: Xi^H_$h($j) = ", ($pre * $E[$n])->bstr, "\n";
