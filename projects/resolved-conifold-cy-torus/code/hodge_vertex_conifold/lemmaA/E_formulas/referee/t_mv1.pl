# t_mv1.pl hmax : the one-partition Marino-Vafa formula (Liu-Liu-Zhou 2003, Thm 1; Okounkov-Pandharipande 2004) as an independent
# test of the Hodge class at every genus.  For mu = (d) the MV formula C_(d) = R_(d) (log of one part = the part) reduces, with
# chi_{(d-k,1^k)}((d)) = (-1)^k, kappa_{hook} = d(d-1-2k), V_hook = 1/([d][d-k-1]![k]!), [m] = 2 sin(m lambda/2), and the q-binomial
# theorem, to (derivation in referee.md):
#   I_g(d) := Int_{Mbar_{g,1}} Lambda^vee(1) Lambda^vee(tau) Lambda^vee(-1-tau) / (1 - d psi)
#           = (1/d^2) [lambda^{2g}]  prod_{a=1}^{d-1} f((d tau + a) lambda) / prod_{a=1}^{d} f(a lambda),   f(y) = sin(y/2)/(y/2),
# (genus 0 term 1/d^2 = the unstable convention).  Here Omega = prod_w Lambda^vee(w)/w^g, so Int Omega lam_d = (-sigma)^{-g} I_g(d).
# Checks: (a) the formula against the Mathematica table mv_g<g>_l1.txt (c-normalisation) for g <= 4 (fixes every normalisation);
#         (b) the referee engine against the formula for g <= hmax, d = 1..6, several tau.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use RefParse; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my $hmax = shift // 5; my $t0 = time;
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; $r }
my $K = $hmax + 1;
sub fser { my $x = shift; return [ map { ((($_ % 2) ? Q(-1) : Q(1)) * $x ** (2*$_)) / (Q(4) ** $_ * fact(2*$_ + 1)) } 0..$K ]; }   # f(x lambda), coeff of lambda^{2k}
sub smul { my ($a, $b) = @_; my @c = map { Q(0) } 0..$K; for my $i (0..$K) { for my $j (0..$K-$i) { $c[$i+$j] += $a->[$i] * $b->[$j]; } } return \@c; }
sub sinv { my $a = shift; my @r = (Q(1) / $a->[0]); for my $n (1..$K) { my $s = Q(0); for my $m (1..$n) { $s += $a->[$m] * $r[$n-$m]; } push @r, -$s / $a->[0]; } return \@r; }
sub formula { my ($g, $d, $tau) = @_; my $num = [Q(1), map { Q(0) } 1..$K]; my $den = [Q(1), map { Q(0) } 1..$K];
  for my $a (1..$d-1) { $num = smul($num, fser($d * $tau + $a)); } for my $a (1..$d) { $den = smul($den, fser(Q($a))); }
  return smul($num, sinv($den))->[$g] / ($d * $d); }
# (a) against the Mathematica tables
my $DATA = "$FindBin::Bin/../../data";
for my $g (1..4) { my $f = "$DATA/mv_g${g}_l1.txt"; next unless -e $f; open my $fh, '<', $f or die; local $/; my $line = <$fh>; close $fh; my ($gg, $ll, $e) = split /\s*\|\s*/, $line, 3;
  my ($n, $bad) = (0, 0);
  for my $x (map { Q($_) } qw(2 1/3 -3/5 5)) { my $c = 1 + $x; my $tau = -$x / $c; my $sig = $tau * (1 + $tau);
    for my $d (1..6) { my $tab = RefParse::evalx($e, { x => $x, m1 => Q($d) }); my $pred = formula($g, $d, $tau) / ((-$sig) ** $g) / ($c ** (3*$g - 1)); $n++; $bad++ if $tab != $pred; } }
  printf "(a) MV one-part formula vs Mathematica table mv_g%d_l1: %d comparisons, %d mismatches\n", $g, $n, $bad; }
# (b) referee engine against the formula
for my $g (1..$hmax) { my ($n, $bad) = (0, 0); my $maxd = $g == 1 ? 1 : 3*$g - 3;
  for my $tau (map { Q($_) } qw(2 -1/3 3/7)) { my $om = RefHodge::omega($g, [Q(1), $tau, -1 - $tau], $maxd); my $sig = $tau * (1 + $tau);
    for my $d (1..6) { my $q = RefHodge::xi($g, $om, [$d], 0); my $mine = $q->[0] * ((-$sig) ** $g); my $pred = formula($g, $d, $tau); $n++; $bad++ if $mine != $pred; } }
  printf "(b) genus %d: referee engine vs MV one-part formula: %d comparisons, %d mismatches [%.1f s]\n", $g, $n, $bad, time - $t0; }
