#!/usr/bin/perl
# r = 1: compare the Toda recursion (Rec.pm) with the closed formula
#   Xi^(0)_g(j) = (j-1)(j-2)...(j-2g+2) [t^{2g}] S(t)^2 H(t)^j = (j-1)...(j-2g+2) [y^j] x(y)^{j-2g},  x = y H(x),
#   S(t) = (t/2)/sinh(t/2),  log H = sum_k (2k-1) B_{2k} t^{2k}/(2k (2k)!).
# Usage: perl closed_r1_check.pl GMAX
use strict; use warnings; use lib '.'; use Rec; use Math::BigRat; use Math::BigInt;
sub BR { Math::BigRat->new(shift) }
my $GMAX = shift // 7;
# Bernoulli numbers B_0..B_{2GMAX}
my @B = (BR(1)); for my $m (1..2*$GMAX) { my $s = BR(0); for my $k (0..$m-1) { $s += BR(binom($m+1, $k)) * $B[$k] } push @B, -$s / BR($m+1) }
sub binom { my ($n, $k) = @_; my $r = Math::BigInt->new(1); for my $i (1..$k) { $r *= ($n - $k + $i); $r /= $i } return $r }
sub fact { my $n = shift; my $r = Math::BigInt->new(1); $r *= $_ for 1..$n; return $r }
my @logS = (BR(0)); my @logH = (BR(0));      # coefficients of t^{2k}
for my $k (1..$GMAX) { my $c = $B[2*$k] / (BR(2*$k) * BR(fact(2*$k))); push @logS, -$c; push @logH, BR(2*$k - 1) * $c }
sub expser { my @a = @_; my @r = (BR(1)); for my $n (1..$#a) { my $s = BR(0); for my $m (1..$n) { $s += BR($m) * $a[$m] * $r[$n - $m] } push @r, $s / BR($n) } return @r }
my ($n, $bad) = (0, 0);
for my $g (1..$GMAX) { for my $j (2..3*$g + 2) {
  my @ex = map { 2 * $logS[$_] + BR($j) * $logH[$_] } 0..$GMAX; my @ser = expser(@ex);
  my $closed = Rec::ff($j - 1, 2*$g - 2) * $ser[$g];
  my $delta = $j - (2*$g - 1);
  if ($delta == 0) { printf "g=%d j=%d (delta=0): closed formula %s (recursion silent)\n", $g, $j, $closed; next }
  my $rec = Rec::xi($g, [$j]); $n++; my $ok = $rec == $closed ? "ok" : "MISMATCH"; $bad++ if $ok ne "ok";
  printf "g=%d j=%d: recursion %s, closed %s %s\n", $g, $j, $rec, $closed, $ok } }
print "closed_r1_check up to genus $GMAX: $n comparisons, $bad mismatches\n";
