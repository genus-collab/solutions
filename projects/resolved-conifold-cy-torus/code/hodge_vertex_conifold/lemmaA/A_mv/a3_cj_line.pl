#!/usr/bin/perl
# a3_cj_line.pl -- check of the cut-and-join identity on the unit-leg line (derived from the Marino-Vafa cut-and-join):
#     (g - 1 + a + sigma d/dsigma) M_g(1^a) = a(a-1) M_g(2, 1^(a-2)),     M_g(mu) = <Omega prod lam_{mu_i}>_g  (c = 1),
# for g >= 2 (a = 0 uses <Omega>_g ~ sigma^(1-g) by the CY degree argument; not computed here: a >= 1).
# sigma^g M_g is a polynomial of degree <= g in sigma; it is interpolated exactly from g+2 values of tau (one extra value as a check).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($g, $amax) = @ARGV; $g //= 2; $amax //= 6;
my @taus = map { R($_) } ('1', '2', '3', '1/2', '5', '7/3', '4', '3/2')[0..$g+1];
sub Mg { my ($mu, $tau) = @_; my $l = scalar @$mu; my $n = 0; $n += $_ for @$mu; my $sig = $tau*(1+$tau);
  my $s = mv_conn($mu, $tau, 2*$g - 2 + $l + $n);
  my $Pi = R(1); for my $j (@$mu) { for my $b (1..$j-1) { $Pi *= ($j*$tau + $b); } $Pi /= fact($j-1); }
  my $I = R((-1)**$g) * $s->[2*$g - 2 + $l + $n] * mvlib::autf(@$mu) / ($sig**($l-1) * $Pi); return $I / ((-$sig)**$g); }
sub interp { my @pts = @_; my $n = scalar @pts; my @c = map { R(0) } 0..$n-1;
  for my $i (0..$n-1) { my @b = (R(1)); my $den = R(1);
    for my $j (0..$n-1) { next if $j == $i; my @nb = map { R(0) } 0..scalar(@b); for my $k (0..$#b) { $nb[$k+1] += $b[$k]; $nb[$k] -= $b[$k] * $pts[$j][0]; } @b = @nb; $den *= ($pts[$i][0] - $pts[$j][0]); }
    for my $k (0..$#b) { $c[$k] += $b[$k] * $pts[$i][1] / $den; } } return @c; }
my $bad = 0;
for my $a (1..$amax) {
  # f(sigma) = sigma^g M_g(1^a) as polynomial; (g-1+a+sigma d/dsigma) M = sigma^{-g} * sum_k (k - 1 + a) c_k sigma^k
  my @p1 = map { my $sig = $_*(1+$_); [ $sig, $sig**$g * Mg([(1) x $a], $_) ] } @taus; my @c1 = interp(@p1);
  my @lhs = map { ($_ - 1 + $a) * $c1[$_] } 0..$#c1;
  my @rhs = map { R(0) } 0..$#c1;
  if ($a >= 2) { my @p2 = map { my $sig = $_*(1+$_); [ $sig, $sig**$g * Mg([2, (1) x ($a-2)], $_) ] } @taus; my @c2 = interp(@p2); @rhs = map { $a*($a-1) * $c2[$_] } 0..$#c2; }
  my $ok = 1; for my $k (0..$#lhs) { $ok = 0 unless $lhs[$k] == $rhs[$k]; }
  $bad++ unless $ok; my $top = $c1[$#c1];
  printf "g=%d a=%d: %s   (sigma^g M_g(1^a) = %s; degree check: coefficient of sigma^%d = %s)\n", $g, $a, ($ok ? "identity holds" : "FAILS"),
    join(' + ', map { "($c1[$_])s^$_" } 0..$g), $g+1, $top; }
print "SUMMARY: ", ($bad ? "$bad failures" : "all identities hold"), "\n";
