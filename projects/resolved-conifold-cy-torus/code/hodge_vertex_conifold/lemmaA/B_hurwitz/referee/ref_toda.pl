#!/usr/bin/perl
# Referee check of the specialised Toda equation (report (T); Okounkov, MRL 7 (2000), eq. (10) at P' = (1,0,0,..), q = 1):
#     d/dp_1 ( E H ) = exp( H(e^b o p) - 2 H(p) + H(e^{-b} o p) ),   E = sum k p_k d/dp_k,
# with H the connected single-Hurwitz series computed by the referee's cut-and-join engine (no characters at all here) and the
# integer E-logarithm; exp(.) computed by the E-recursion  |nu| Y(nu) = sum_{m1 | nu} |m1| X(m1) Y(nu/m1).
# Usage: perl ref_toda.pl DMAX BMAX
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefHur; use Math::BigRat; use Math::BigInt;
sub BR { Math::BigRat->new(shift) } sub BI { Math::BigInt->new(shift) }
my ($DMAX, $BMAX) = @ARGV; $DMAX //= 9; $BMAX //= 16;
# all monomials of degree <= DMAX as partitions
my @mons; for my $d (0..$DMAX) { push @mons, @{RefHur::parts_of($d)} }
my %deg; $deg{RefHur::pkey($_)} = do { my $s = 0; $s += $_ for @$_; $s } for @mons;
sub submons { my $mu = shift; my %m; $m{$_}++ for @$mu; my @v = sort { $b <=> $a } keys %m; my @out = ([]);
  for my $x (@v) { my @n; for my $s (@out) { push @n, [@$s, ($x) x $_] for 0..$m{$x} } @out = @n } return @out }
sub mdiv { my ($mu, $m1) = @_; my %m; $m{$_}++ for @$mu; $m{$_}-- for @$m1; return [ sort { $b <=> $a } map { ($_) x $m{$_} } grep { $m{$_} > 0 } keys %m ] }
# Zhat and Hhat (integers) for all monomials
my (%Z, %H);
for my $d (0..$DMAX) { my $T = RefHur::cj_table($d, $BMAX);
  for my $mu (@{RefHur::parts_of($d)}) { my $k = RefHur::pkey($mu); my $cs = RefHur::classsize($mu); $Z{$k} = [ map { $T->{$k}[$_] * $cs } 0..$BMAX ] } }
for my $mu (@mons) { my $k = RefHur::pkey($mu); my $d = $deg{$k}; next if $d == 0; my @h = map { $_->copy } @{$Z{$k}};
  for my $m1 (submons($mu)) { my $k1 = RefHur::pkey([sort { $b <=> $a } @$m1]); my $d1 = 0; $d1 += $_ for @$m1; next if $d1 == 0 || $d1 == $d;
    my $k2 = RefHur::pkey(mdiv($mu, $m1)); my $cd = RefHur::bin($d - 1, $d1 - 1);
    for my $b1 (0..$BMAX) { next if $H{$k1}[$b1]->is_zero; for my $b ($b1..$BMAX) { next if $Z{$k2}[$b - $b1]->is_zero; $h[$b] -= $cd * RefHur::bin($b, $b1) * $H{$k1}[$b1] * $Z{$k2}[$b - $b1] } } }
  $H{$k} = \@h }
# H as BigRat ordinary coefficients: [p_mu beta^b] H = Hhat/(d! b!)
my %Hq; for my $k (keys %H) { my $d = $deg{$k}; $Hq{$k} = [ map { BR($H{$k}[$_]) / (BR(RefHur::fac($d)) * BR(RefHur::fac($_))) } 0..$BMAX ] }
# X = Delta H: multiply [p_mu beta^b] by (2cosh(d beta) - 2)
my %X; for my $k (keys %Hq) { my $d = $deg{$k}; my @x = map { BR(0) } 0..$BMAX;
  for my $b (0..$BMAX) { next if $Hq{$k}[$b]->is_zero; for (my $e = 2; $b + $e <= $BMAX; $e += 2) { $x[$b + $e] += $Hq{$k}[$b] * BR(2) * BR(BI($d) ** $e) / BR(RefHur::fac($e)) } }
  $X{$k} = \@x }
# Y = exp(X) on monomials of degree <= DMAX-1
my %Y = ("" => [BR(1), map { BR(0) } 1..$BMAX]);
for my $nu (@mons) { my $k = RefHur::pkey($nu); my $d = $deg{$k}; next if $d == 0 || $d > $DMAX - 1; my @y = map { BR(0) } 0..$BMAX;
  for my $m1 (submons($nu)) { my $d1 = 0; $d1 += $_ for @$m1; next if $d1 == 0; my $k1 = RefHur::pkey([sort { $b <=> $a } @$m1]); my $k2 = RefHur::pkey(mdiv($nu, $m1));
    for my $b1 (0..$BMAX) { next if $X{$k1}[$b1]->is_zero; for my $b ($b1..$BMAX) { next if $Y{$k2}[$b - $b1]->is_zero; $y[$b] += BR($d1) * $X{$k1}[$b1] * $Y{$k2}[$b - $b1] } } }
  $Y{$k} = [ map { $_ / BR($d) } @y ] }
# LHS: d/dp_1 (E H) at nu: mu = nu u {1}: |mu| m_1(mu) [p_mu] H
my ($n, $bad, $nz) = (0, 0, 0);
for my $nu (@mons) { my $k = RefHur::pkey($nu); my $d = $deg{$k}; next if $d > $DMAX - 1;
  my $mu = [sort { $b <=> $a } (@$nu, 1)]; my $km = RefHur::pkey($mu); my $m1 = grep { $_ == 1 } @$mu;
  for my $b (0..$BMAX) { my $lhs = $Hq{$km}[$b] * BR(($d + 1) * $m1); my $rhs = $Y{$k}[$b];
    $n++; $nz++ unless $lhs->is_zero && $rhs->is_zero; if ($lhs != $rhs) { $bad++; print "MISMATCH p_[$k] beta^$b: LHS $lhs RHS $rhs\n" } } }
print "Toda (T), referee engine: monomials of degree <= ", $DMAX - 1, ", beta^<= $BMAX: $n coefficients ($nz non-zero), $bad mismatches\n";
# the opposite sign of Delta (exp(-Delta H)) must fail: sanity check that the test has teeth
my $fail = 0; for my $nu (@mons) { my $k = RefHur::pkey($nu); next unless $deg{$k} == 1; my $mu = [1, 1];
  my $lhs = $Hq{"1,1"}[4] * BR(2 * 2); my $rhs_minus = BR(0); # [p_1 beta^4] of exp(-X): -X(p_1)[4] (only one term at degree 1)
  $rhs_minus = -$X{"1"}[4]; $fail = 1 if $lhs != $rhs_minus }
print "sanity: sign-flipped right side fails at p_1 beta^4: ", ($fail ? "yes" : "NO"), "\n";
