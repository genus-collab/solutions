#!/usr/bin/perl
# ref_limits.pl HMIN HMAX : the two end components of (6.1) against the independent all-genus TR-free theorems:
#  (a) Hurwitz component (k = 0): Route B, Theorem B (ELSV + Toda):  Xi^(0)_h(j) = (j-1)...(j-2h+2) [t^{2h}] S(t)^2 H(t)^j,  H = S exp(-t S'/S);
#  (b) lambda_g component (k = h):  sigma^0-coefficient of sigma^h Xi_h(j) = (-1)^h Xi^(h)_h(j) must equal prod_{l=1}^{2h-2}(j-l) [t^{2h}] S(t)
#      (Faber-Pandharipande lambda_g formula in the normalisation of brief s.2(a)).
# Both are checked as exact polynomial identities in j.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($h1, $h2) = @ARGV;
# univariate series in u = t^2 with polynomial-in-j coefficients: array of hashes b => coef (j^b)
sub jmul { my ($p, $q) = @_; my %r; for my $a (keys %$p) { for my $b (keys %$q) { my $k = $a + $b; my $x = $p->{$a} * $q->{$b}; if (exists $r{$k}) { $r{$k} += $x } else { $r{$k} = $x } } } return \%r; }
sub jadd { my ($p, $q, $c) = @_; $c //= Q(1); my %r = map { $_ => $p->{$_}->copy } keys %$p; for (keys %$q) { my $x = $q->{$_} * $c; if (exists $r{$_}) { $r{$_} += $x } else { $r{$_} = $x } } return \%r; }
for my $h ($h1..$h2) {
  my @c = (Q(0)); for my $k (1..$h) { push @c, -RefR2::bern(2*$k) / (Q(2*$k) * RefR2::fact(2*$k)); }   # log S = Sum c_k t^{2k}
  # exponent of S^2 H^j = exp( 2 log S + j (log S - t (log S)') ) = exp( Sum_k c_k u^k (2 + j (1 - 2k)) )
  my @e = ({}); for my $k (1..$h) { push @e, { 0 => 2 * $c[$k], 1 => (1 - 2*$k) * $c[$k] }; }
  my @f = ({ 0 => Q(1) }); for my $n (1..$h) { my $acc = {}; for my $k (1..$n) { $acc = jadd($acc, jmul($e[$k], $f[$n-$k]), Q($k)); } my %t = map { $_ => $acc->{$_} / $n } keys %$acc; push @f, \%t; }
  my $P = $f[$h]; for my $l (1..2*$h-2) { $P = jmul($P, { 1 => Q(1), 0 => Q(-$l) }); }
  my $C = RefR2::comps61($h); my $bad = 0; my %keys = map { $_ => 1 } (keys %$P, keys %{$C->{0}});
  for my $b (keys %keys) { my $x = $P->{$b} // Q(0); my $y = $C->{0}{$b} // Q(0); $bad++ if $x != $y; }
  # lambda_g: b_h = [t^{2h}] (t/2)/sin(t/2) by series inversion of sin(t/2)/(t/2) (no Bernoulli numbers); expected sigma^0 part = (-1)^h b_h prod(j-l)
  my @sn = map { (($_ % 2) ? Q(-1) : Q(1)) / (Q(4)**$_ * RefR2::fact(2*$_ + 1)) } 0..$h; my @g = (Q(1)); for my $n (1..$h) { my $s = Q(0); $s -= $sn[$_] * $g[$n-$_] for 1..$n; push @g, $s; }
  my $L = { 0 => $g[$h] * (($h % 2) ? Q(-1) : Q(1)) }; for my $l (1..2*$h-2) { $L = jmul($L, { 1 => Q(1), 0 => Q(-$l) }); }
  my $bad2 = 0; my $sg = ($h % 2) ? Q(-1) : Q(1); my %k2 = map { $_ => 1 } (keys %$L, keys %{$C->{$h}});
  for my $b (keys %k2) { my $x = $L->{$b} // Q(0); my $y = $sg * ($C->{$h}{$b} // Q(0)); $bad2++ if $x != $y; }
  printf "h=%d: Hurwitz component of (6.1) vs Route B Theorem B: %d coefficient mismatches; lambda_g component vs lambda_g formula: %d mismatches\n", $h, $bad, $bad2; }
