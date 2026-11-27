# check61_engine.pl h1 h2 ... : compare the one-leaf closed formula (6.1) (proof_M_residue.md s.6)
#   sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j - l) [lambda^{2h}] S(e1 l)^{1 + j e2} S(e2 l)^{1 + j e1},  e1 = -tau, e2 = 1 + tau
# with the engine components (engine_h<h>_r1.txt), identically in j and tau (by comparing all sigma-coefficients through
# enough evaluation points: deg_j <= 3h-2, deg_{1/sigma} <= h).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Comps;
sub R { Math::BigRat->new(shift) }
my @B = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"));
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
sub formula { my ($tau, $j, $h) = @_; my $e1 = -$tau; my $e2 = 1 + $tau; my $K = $h;
  my @c = (R(0)); for my $k (1..$K) { my $lS = -$B[$k] / (2*$k * fact(2*$k)); push @c, (1 + $j*$e2) * $lS * $e1**(2*$k) + (1 + $j*$e1) * $lS * $e2**(2*$k); }
  my @r = (R(1)); for my $n (1..$K) { my $s = R(0); for my $m (1..$n) { $s += $m * $c[$m] * $r[$n-$m] } $r[$n] = $s / $n; }
  my $pre = R(1); $pre *= ($j - $_) for 1..(2*$h-2); my $sig = $tau * (1 + $tau);
  return $pre * $r[$h] / $sig**$h; }
for my $h (@ARGV) {
  my $C = Comps::load("$FindBin::Bin/engine_h${h}_r1.txt");
  my ($n, $bad) = (0, 0);
  for my $tau (map { R($_) } ("1/3", "2", "-3/7", "5", "-1/4", "7/2", "1/9")) {
    for my $j (0..(3*$h)) { $n++; my $a = Comps::xi_at($C, $h, 1, [$j], $tau); my $b = formula($tau, R($j), $h); $bad++ if $a != $b; } }
  print "genus $h: (6.1) vs engine at $n points (7 values of tau, j = 0..3h): $bad mismatches\n"; }
