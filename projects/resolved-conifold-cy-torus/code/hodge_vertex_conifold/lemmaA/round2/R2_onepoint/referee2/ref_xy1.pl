#!/usr/bin/perl
# ref_xy1.pl TAU GMAX JMAX : referee-2 implementation (from scratch) of Theorem 1 at n = 1 for an integer framing TAU >= 1,
#   [x^j] F_{g,1} = j^{-2} [hb^{2g}] [v^{-1}]_{v=oo} Q_j(v),
#   Q_j(v) = prod_{m in M(e2 j)} (v + 1 + m hb/e2) / prod_{m in M(tau j)} (v + m hb/tau),   M(N) = {k - (N+1)/2 : k = 1..N},
# then Xi_g(j) = sigma^{g-1} [x^j] F_{g,1} / Pi_j (dictionary P2), compared EXACTLY with formula (6.1) (RefR2::rhs61).
# Implementation: w = 1/v, [v^{-1}] Q_j = [w^{j+1}] prod (1 + (1 + m hb/e2) w) / prod (1 + (m hb/tau) w), a truncated series in (w, hb).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($tau, $GMAX, $JMAX) = @ARGV; $tau = Q($tau); my $e2 = $tau + 1; my $sig = $tau * $e2;
die "integer tau >= 1 only" unless $tau->is_int && $tau >= 1;
my %cache; my ($n, $bad, $nz) = (0, 0, 0);
for my $j (1..$JMAX) {
  my $W = $j + 1; my $HB = 2 * $GMAX;
  # F[w][hb]
  my @F = map { [ map { Q(0) } 0..$HB ] } 0..$W; $F[0][0] = Q(1);
  my $N2 = $e2 * $j; my $N1 = $tau * $j;
  for my $kk (1..$N2) { my $m = Q($kk) - ($N2 + 1) / 2; my $c1 = $m / $e2;   # factor 1 + w + c1 hb w
    for (my $a = $W; $a >= 1; $a--) { for my $p (0..$HB) { my $add = $F[$a-1][$p]; $add += $c1 * $F[$a-1][$p-1] if $p >= 1; $F[$a][$p] += $add; } } }
  for my $kk (1..$N1) { my $m = Q($kk) - ($N1 + 1) / 2; my $c = $m / $tau;    # divide by 1 + c hb w : G_a = F_a - c hb G_{a-1}
    for my $a (1..$W) { for (my $p = $HB; $p >= 1; $p--) { $F[$a][$p] -= $c * $F[$a-1][$p-1]; } } }
  for my $g (1..$GMAX) {
    my $Fx = $F[$W][2*$g] / ($j * $j);
    my $xi = $sig**($g-1) * $Fx / RefR2::Pi($j, $tau);
    my $pred = RefR2::xi61($g, $tau, $j, \%cache);
    $n++; $nz++ unless $xi->is_zero; if ($xi != $pred) { $bad++; print "MISMATCH tau=$tau g=$g j=$j: xy $xi  (6.1) $pred\n"; }
  } }
print "ref_xy1 tau=$tau g<=$GMAX j<=$JMAX: compared $n, non-zero $nz, mismatches $bad\n";
