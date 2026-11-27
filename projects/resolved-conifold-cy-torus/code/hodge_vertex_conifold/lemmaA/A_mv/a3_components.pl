#!/usr/bin/perl
# a3_components.pl -- tau-independent components of the Marino-Vafa integrals:
#   I_g(mu; tau) = Int Lambda(1)Lambda(tau)Lambda(-1-tau) prod 1/(1 - mu_i psi_i) = Sum_{k=0}^{g} (-1)^{g-k} sigma^{g-k} m_k(g, mu),
#   m_k(g, mu) = Int Lambda(1) e_k(alpha + alpha^2) prod 1/(1 - mu_i psi_i)   (k = g: lambda_g formula, k = 0: Hurwitz/ELSV).
# Computed from the Schur side at g+1 values of tau (exact Lagrange interpolation in sigma); prints m_k for the requested mu.
# Usage: perl a3_components.pl g "mu1;mu2;..."      (mu as comma lists)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($g, $mulist) = @ARGV; my @mus = map { [ split /,/, $_ ] } split /;/, $mulist;
my @taus = map { R($_) } ('1', '2', '3', '1/2', '5', '7/3', '4');
@taus = @taus[0..$g];   # g+1 values, distinct sigma
sub Ivals { my ($mu) = @_; my $l = scalar @$mu; my $n = 0; $n += $_ for @$mu; my @v;
  for my $tau (@taus) { my $sig = $tau*(1+$tau); my $P = 2*$g - 2 + $l + $n; my $s = mv_conn($mu, $tau, $P);
    my $Pi = R(1); for my $j (@$mu) { for my $b (1..$j-1) { $Pi *= ($j*$tau + $b); } $Pi /= fact($j-1); }
    push @v, [ $sig, R((-1)**$g) * $s->[2*$g - 2 + $l + $n] * mvlib::autf(@$mu) / ($sig**($l-1) * $Pi) ]; }
  return @v; }
# interpolate polynomial in sigma of degree <= g through points
sub interp { my @pts = @_; my $n = scalar @pts; my @c = map { R(0) } 0..$n-1;
  for my $i (0..$n-1) { my @b = (R(1)); my $den = R(1);
    for my $j (0..$n-1) { next if $j == $i; my @nb = map { R(0) } 0..scalar(@b); for my $k (0..$#b) { $nb[$k+1] += $b[$k]; $nb[$k] -= $b[$k] * $pts[$j][0]; } @b = @nb; $den *= ($pts[$i][0] - $pts[$j][0]); }
    for my $k (0..$#b) { $c[$k] += $b[$k] * $pts[$i][1] / $den; } } return @c; }
for my $mu (@mus) { my @c = interp(Ivals($mu));   # I_g = sum_j c_j sigma^j ; c_j = (-1)^j m_{g-j}
  my @m = map { R((-1)**($g - $_)) * $c[$g - $_] } 0..$g;  # m[k] = (-1)^{g-k} c_{g-k}
  printf "g=%d mu=(%s): %s\n", $g, join(',', @$mu), join('  ', map { "m$_=" . $m[$_] } 0..$g); }
