#!/usr/bin/perl
# Referee check of the specialised Toda equation (report (T); Okounkov, MRL 7 (2000), eq. (10) at P' = (1,0,0,..), q = 1):
#     d/dp_1 ( E H ) = exp( H(e^b o p) - 2 H(p) + H(e^{-b} o p) ),   E = sum k p_k d/dp_k,
# with H the connected single-Hurwitz series computed by the referee's cut-and-join engine (no characters at all here) and the
# integer E-logarithm; exp(.) computed by the E-recursion  |nu| Y(nu) = sum_{m1 | nu} |m1| X(m1) Y(nu/m1).
# Usage: perl ref_mkp.pl DMAX BMAX
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
# ---- referee check of the report's modified-KP relation (claims list, computer-verified for the Hurwitz series):
#   s * 2 d_{p_2}(F - F1) + d_{p_1}^2 (F + F1) + (d_{p_1}(F - F1))^2 = 0,  F = H(p), F1 = H(e^beta o p), s = +1 (s = -1 must fail)
my %F1; for my $k (keys %Hq) { my $d = $deg{$k}; my @x = map { BR(0) } 0..$BMAX;
  for my $b (0..$BMAX) { next if $Hq{$k}[$b]->is_zero; for my $e (0..$BMAX - $b) { $x[$b + $e] += $Hq{$k}[$b] * BR(BI($d) ** $e) / BR(RefHur::fac($e)) } } $F1{$k} = \@x }
sub comb { my ($A, $B, $c) = @_; my %S; for my $k (keys %$A) { $S{$k} = [ map { $A->{$k}[$_] + $c * $B->{$k}[$_] } 0..$BMAX ] } \%S }
sub dp { my ($A, $j) = @_; my %S;   # d/dp_j on monomials (partitions)
  for my $k (keys %$A) { my @mu = $k eq '' ? () : split /,/, $k; my $m = grep { $_ == $j } @mu; next unless $m;
    my @rest = @mu; for my $i (0..$#rest) { if ($rest[$i] == $j) { splice(@rest, $i, 1); last } } my $kr = join(",", @rest);
    $S{$kr} //= [ map { BR(0) } 0..$BMAX ]; $S{$kr}[$_] += BR($m) * $A->{$k}[$_] for 0..$BMAX } \%S }
sub smul2 { my ($A, $B, $dmax) = @_; my %S;
  for my $ka (keys %$A) { for my $kb (keys %$B) { my @m = sort { $b <=> $a } ((split /,/, $ka), (split /,/, $kb)); my $d = 0; $d += $_ for @m; next if $d > $dmax;
      my $k = join(",", @m); $S{$k} //= [ map { BR(0) } 0..$BMAX ];
      for my $i (0..$BMAX) { next if $A->{$ka}[$i]->is_zero; for my $j (0..$BMAX - $i) { $S{$k}[$i + $j] += $A->{$ka}[$i] * $B->{$kb}[$j] } } } } \%S }
my $Dm = comb(\%Hq, \%F1, BR(-1)); my $Sp = comb(\%Hq, \%F1, BR(1));
my $T2 = dp($Dm, 2); my $T11 = dp(dp($Sp, 1), 1); my $T1 = dp($Dm, 1); my $Q = smul2($T1, $T1, $DMAX - 2);
for my $s (1, -1) { my ($n, $nz, $bad) = (0, 0, 0); my %keys; $keys{$_} = 1 for (keys %$T2, keys %$T11, keys %$Q);
  for my $k (keys %keys) { my $d = 0; $d += $_ for ($k eq '' ? () : split /,/, $k); next if $d > $DMAX - 2;
    for my $b (0..$BMAX) { my $v = BR(0); $v += $T11->{$k}[$b] if $T11->{$k}; $v += $Q->{$k}[$b] if $Q->{$k}; $v += BR(2 * $s) * $T2->{$k}[$b] if $T2->{$k};
      $n++; $bad++ unless $v->is_zero } }
  print "mKP relation, sign s = $s: $n coefficients (degree <= ", $DMAX - 2, ", beta^<= $BMAX), non-zero: $bad\n" }
