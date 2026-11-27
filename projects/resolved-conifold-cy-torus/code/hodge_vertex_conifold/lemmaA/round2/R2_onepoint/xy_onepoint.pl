#!/usr/bin/perl
# xy_onepoint.pl tau jmax G : compare, for integer tau >= 1, the x-y swap closed formula
#   (*)  sigma^g Xi_g(j) = (1/(sigma j^2 Pi_j)) [lambda^{2g}] [z^{-1}]_{z=oo}
#            prod_{m' in M2} (z + 1 + m' tau lambda) / prod_{m in M1} (z + m e2 lambda),
#        M1 = {k - (tau j + 1)/2 : k = 1..tau j},  M2 = {k - (e2 j + 1)/2 : k = 1..e2 j},  e2 = 1 + tau, sigma = tau e2,
# with the one-leaf formula (6.1):
#   (6.1) sigma^g Xi_g(j) = prod_{l=1}^{2g-2} (j - l) [lambda^{2g}] S(tau lambda)^{1 + j e2} S(e2 lambda)^{1 - j tau}.
# Exact arithmetic (Math::BigRat).  [z^{-1}] at z = oo:  with t = 1/z,  Q = z^j * prod(1 + (1 + m' tau l) t) / prod(1 + m e2 l t),
# so [z^{-1}] Q = [t^{j+1}] of  prod(1 + (1 + m' tau l) t) * prod(1 + m e2 l t)^{-1}.
use strict; use warnings;
use Math::BigRat;
sub Q { Math::BigRat->new(@_) }
my ($tau, $jmax, $G) = @ARGV; $tau //= 1; $jmax //= 8; $G //= 4;
$tau = Q($tau); my $e2 = $tau + 1; my $sig = $tau * $e2;
my $L = 2*$G;   # lambda truncation
# bivariate series: array over t-degree (0..T) of arrays over lambda-degree (0..L)
sub bz { my ($T) = @_; return [ map { [ map { Q(0) } 0..$L ] } 0..$T ]; }
sub bmul { my ($A, $B, $T) = @_; my $C = bz($T);
  for my $i (0..$T) { for my $k (0..$T-$i) { my $a = $A->[$i]; my $b = $B->[$k]; my $c = $C->[$i+$k];
      for my $p (0..$L) { next if $a->[$p]->is_zero; for my $q (0..$L-$p) { next if $b->[$q]->is_zero; $c->[$p+$q] += $a->[$p] * $b->[$q]; } } } }
  return $C; }
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; return $r; }
sub bern { my $n = shift; my @B = (Q(1)); for my $m (1..$n) { my $s = Q(0); for my $k (0..$m-1) { $s += binomq($m+1, $k) * $B[$k]; } $B[$m] = -$s / ($m+1); } return $B[$n]; }
sub binomq { my ($n, $k) = @_; return Q(0) if $k < 0 || $k > $n; return fact($n) / (fact($k) * fact($n-$k)); }
# log S(y) = Sum_n c_n y^{2n}, c_n = -B_{2n}/(2n (2n)!)
my @c = (Q(0)); for my $n (1..$G) { push @c, -bern(2*$n) / (2*$n * fact(2*$n)); }
sub series_exp { my ($e) = @_; my @f = (Q(1)); for my $n (1..$L) { my $s = Q(0); for my $m (1..$n) { $s += $m * $e->[$m] * $f[$n-$m]; } $f[$n] = $s / $n; } return \@f; }
my ($nbad, $ntot) = (0, 0);
for my $j (1..$jmax) {
  my $T = $j + 1;
  my $n1 = $tau * $j; my $n2 = $e2 * $j;
  # product over M2 of (1 + (1 + m' tau l) t)
  my $P = bz($T); $P->[0][0] = Q(1);
  for my $k (1..$n2) { my $m = Q($k) - ($n2 + 1) / 2; my $F = bz($T); $F->[0][0] = Q(1); $F->[1][0] = Q(1); $F->[1][1] = $m * $tau; $P = bmul($P, $F, $T); }
  # product over M1 of (1 + m e2 l t)^{-1} = sum_k (-m e2 l t)^k
  for my $k (1..$n1) { my $m = Q($k) - ($n1 + 1) / 2; my $F = bz($T); my $cf = Q(1);
    for my $d (0..$T) { $F->[$d][$d] = $cf if $d <= $L; $cf *= -$m * $e2; } $P = bmul($P, $F, $T); }
  my $Tser = $P->[$T];   # [t^{j+1}], series in lambda
  my $Pi = Q(1); for my $b (1..$j-1) { $Pi *= ($j * $tau + $b); } $Pi /= fact($j-1);
  # (6.1): exponent series E = (1 + j e2) log S(tau l) + (1 - j tau) log S(e2 l)
  my @E = map { Q(0) } 0..$L;
  for my $n (1..$G) { $E[2*$n] = $c[$n] * ((1 + $j*$e2) * $tau**(2*$n) + (1 - $j*$tau) * $e2**(2*$n)); }
  my $F61 = series_exp(\@E);
  for my $g (0..$G) {
    my $lhs = $Tser->[2*$g] / ($sig * $j * $j * $Pi);
    my $rhs;
    if ($g == 0) { $rhs = Q(1) / ($j * ($j + 1)); }
    else { my $pre = Q(1); $pre *= ($j - $_) for 1..(2*$g-2); $rhs = $pre * $F61->[2*$g]; }
    $ntot++; my $ok = ($lhs == $rhs); $nbad++ unless $ok;
    printf "tau=%s j=%d g=%d  xy: %s   (6.1): %s  %s\n", $tau, $j, $g, $lhs, $rhs, $ok ? "OK" : "MISMATCH";
  } }
print "total $ntot, mismatches $nbad\n";
