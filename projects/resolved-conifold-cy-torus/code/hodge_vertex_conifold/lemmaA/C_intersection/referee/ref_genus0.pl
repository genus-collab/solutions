#!/usr/bin/perl
# Referee check of Corollary H0: on Mbar_{0,n}, Xi_0(R) = Sum_m 1/m! < prod_i lam_{j_i} prod_k T(psi_k) >_0 = (J-1)(J-2)...(J-n+3).
# Direct evaluation with the genus-0 numbers <tau_{b_1} ... tau_{b_N}>_0 = (N-3)!/prod b_i!  (no recursion of any kind):
# for leaves with Borel transforms e^{j t} (lam_j) and -(e^t - 1 - t) (T), the m-th term is (n+m-3)!/m! [t^{n+m-3}] e^{Jt} (1+t-e^t)^m.
# The value depends on R only through n and J; tested for n = 3..10 and J = -4..16 (all integers, so also J <= 0 and D <= -3 or not).
use strict; use warnings; use Math::BigRat; use Math::BigInt;
sub fct { Math::BigInt->new(shift)->bfac() }
my $bad = 0; my $tot = 0;
for my $n (3 .. 10) { my $N = $n - 3;
  for my $J (-4 .. 16) {
    my $sum = Math::BigRat->new(0);
    for my $m (0 .. $N) { # T has degree >= 2: 2m <= N + m
      my $deg = $N + $m;
      # series e^{Jt} and (1 + t - e^t)^m to order t^deg
      my @e = map { Math::BigRat->new($J) ** $_ / Math::BigRat->new(fct($_)->bstr) } 0 .. $deg;
      my @g = (Math::BigRat->new(0), Math::BigRat->new(0), map { Math::BigRat->new(-1) / Math::BigRat->new(fct($_)->bstr) } 2 .. $deg); # 1 + t - e^t
      my @p = (Math::BigRat->new(1), map { Math::BigRat->new(0) } 1 .. $deg);
      for (1 .. $m) { my @q = map { Math::BigRat->new(0) } 0 .. $deg; for my $a (0 .. $deg) { next if $p[$a]->is_zero; for my $b (0 .. $deg - $a) { $q[$a + $b] += $p[$a] * $g[$b]; } } @p = @q; }
      my $c = Math::BigRat->new(0); for my $a (0 .. $deg) { $c += $e[$a] * $p[$deg - $a]; }
      $sum += $c * Math::BigRat->new(fct($deg)->bstr) / Math::BigRat->new(fct($m)->bstr); }
    my $pred = Math::BigRat->new(1); $pred *= Math::BigRat->new($J - $_) for 1 .. $N;
    $tot++; $bad++ unless $sum == $pred;
    print "n=$n J=$J: direct = ", $sum->bstr, "  formula = ", $pred->bstr, ($sum == $pred ? '' : '  MISMATCH'), "\n" if ($n <= 5 || $sum != $pred); } }
print "SUMMARY: $tot (n, J) pairs, $bad mismatches\n";
