#!/usr/bin/perl
# t_genus0_curve.pl -- genus-0 test of the "deformed mirror curve" picture behind the closed form of the polar-order relations.
# Genus-0 inverse Lax function at a point P of the slice (P_1 = 0, string equation for unit parts):
#   Y_P(z) = z - sum_{n>=1} z^{-n} sum_T Y_0(T+{1,n}) P^T/|Aut T|,   Y_0(S) = sigma Pi_S Gamma(J_S)/Gamma(J_S - r_S + 3).
# Evaluate along the leaf curve P_j(t) = (-1)^j (Pi_j/j) t^{j-1} rho_j(t) (j >= 2), invert lambda + tau = Y_P(z) for u_P(lambda) = 1/z,
# and compare with the prediction  u_t(lambda) = (Lam - c_t)^tau Lam^{-1-tau} (1 + t Lam)^{-tau} e^{-sigma t},  Lam = (e^{t lambda}-1)/t, c_t = (e^t-1)/t.
# Everything as a double series: orders t^0..t^M, and in 1/lambda (u_P = u(lambda) * (1 + sum_m t^m c_m(lambda))).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M, $NL) = @ARGV; $taus //= '2/7'; $M //= 3; $NL //= 10; my $tau = R($taus); my $sig = $tau*(1+$tau);
sub Pi { Pi_at($_[0], $tau) }
sub X0 { my @S = @_; my $J = 0; $J += $_ for @S; my $r = @S; return fact($J-1) / fact($J - $r + 2); }
# leaf curve coefficients P_j(t) as series in t to order M (j <= M+1 contributes up to t^M)
my $Tt = T_of_t_series($M+1); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $M+1), $Tt, $M+1);
my %Pj;
for my $j (2..$M+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $M+1), $Tt, $M+1);
  my $p = [ R(1), map { R(0) } 1..$M+1 ]; $p = series_mul($p, $h1, $M+1) for 1..$j; my $rho = series_mul($hj, series_inv($p, $M+1), $M+1);
  my $c = R((-1)**$j) * Pi($j) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $c * $rho->[$_-$j+1] : R(0) } 0..$M ]; }
# G(n) := sum_T Y_0(T+{1,n}) P(t)^T/|Aut T| as a series in t (T multisets of parts >= 2 with level <= M)
sub multisets { my ($maxlev) = @_; my @out = ([]); for my $L (1..$maxlev) { for my $q (parts($L, $L)) { push @out, [ map { $_ + 1 } @$q ]; } } return @out; }
my @Ts = multisets($M);
my %Gn;
for my $n (1..$NL) { my @g = map { R(0) } 0..$M;
  for my $T (@Ts) { my $pr = [ R(1), map { R(0) } 1..$M ]; for my $t (@$T) { $pr = series_mul($pr, $Pj{$t}, $M); }
    my $Y0 = $sig; $Y0 = $Y0 * Pi($_) for (@$T, $n); $Y0 = $Y0 * X0(@$T, 1, $n) / aut(@$T);
    for my $i (0..$M) { $g[$i] = $g[$i] + $Y0 * $pr->[$i]; } }
  $Gn{$n} = \@g; }
# Y_P(z) = z - sum_n z^{-n} G(n).  Work with w = 1/z: lambda + tau = 1/w - sum_n w^n G(n)  -> solve for w(lambda) as series in mu = 1/lambda
# Represent functions of mu: arrays over t-order of arrays over mu-power (0..NL).
# Solve iteratively: w = 1/(lambda + tau + sum_n w^n G(n)) = mu / (1 + mu (tau + sum_n w^n G(n)))
my $NM = $NL;
sub bimul { my ($a, $b) = @_; my @c = map { [ map { R(0) } 0..$NM ] } 0..$M; for my $i (0..$M) { for my $j (0..$M-$i) { for my $p (0..$NM) { next if $a->[$i][$p]->is_zero; for my $q (0..$NM-$p) { $c[$i+$j][$p+$q] = $c[$i+$j][$p+$q] + $a->[$i][$p]*$b->[$j][$q]; } } } } return \@c; }
my $zero = sub { [ map { [ map { R(0) } 0..$NM ] } 0..$M ] };
my $w = $zero->(); $w->[0][1] = R(1);   # initial mu
for my $it (1..$NM+2) {
  my $S = $zero->(); $S->[0][0] = $tau;
  my $wp = $zero->(); $wp->[0][0] = R(1);
  for my $n (1..$NL) { $wp = bimul($wp, $w); for my $i (0..$M) { for my $j (0..$M-$i) { next if $Gn{$n}[$j]->is_zero; for my $p (0..$NM) { $S->[$i+$j][$p] = $S->[$i+$j][$p] + $Gn{$n}[$j] * $wp->[$i][$p]; } } } }
  # denominator 1 + mu S
  my $den = $zero->(); $den->[0][0] = R(1); for my $i (0..$M) { for my $p (0..$NM-1) { $den->[$i][$p+1] = $den->[$i][$p+1] + $S->[$i][$p]; } }
  # invert den as a bivariate series (constant term 1)
  my $inv = $zero->(); $inv->[0][0] = R(1); my $x = $zero->(); for my $i (0..$M) { for my $p (0..$NM) { $x->[$i][$p] = -$den->[$i][$p]; } } $x->[0][0] = R(0);
  my $pw = $zero->(); $pw->[0][0] = R(1); for my $k (1..$NM + $M) { $pw = bimul($pw, $x); for my $i (0..$M) { for my $p (0..$NM) { $inv->[$i][$p] = $inv->[$i][$p] + $pw->[$i][$p]; } } }
  my $muS = $zero->(); $muS->[0][1] = R(1); $w = bimul($muS, $inv); }
# prediction: u_t(lambda) as a bivariate series in (t, mu): compute log(u_t/u) = tau log((Lam - c_t)/(lambda-1)) - (1+tau) log(Lam/lambda) - tau log(1 + t Lam) - sigma t
# with Lam = (e^{t lambda} - 1)/t = lambda (1 + t lambda/2 + ...): each t-order is a polynomial in lambda -> negative powers of mu. Compare instead
# the t-coefficients of log(u_P/u) at fixed lambda numerically at several lambda values (large), as truncated Laurent series in mu.
# Here: evaluate both at lambda = Lval (rational) using the mu-series truncated at order NM (u_P) and exact formulas (prediction).
my $Lval = R(shift // 50); my $muv = R(1)/$Lval;
my @uP = map { my $i = $_; my $s = R(0); my $pw = R(1); for my $p (0..$NM) { $s = $s + $w->[$i][$p]*$pw; $pw = $pw*$muv; } $s } 0..$M;   # u_P(lambda) t-coefficients
# u(lambda) = (lambda-1)^tau lambda^{-1-tau} = mu (1 - mu)^tau
my $u0 = $muv * exp_rat_pow(R(1) - $muv, $tau, 40);
sub exp_rat_pow { my ($x, $e, $N) = @_; # (x)^e for x near 1 via binomial series
  my $d = $x - 1; my $s = R(0); my $term = R(1); for my $k (0..$N) { $s = $s + $term; $term = $term * ($e - $k) / ($k + 1) * $d; } return $s; }
print "u_P t-coefficients at lambda = $Lval: ", join(', ', map { sprintf("%.12e", $_->numify) } @uP), "\n";
print "u(lambda) = ", sprintf("%.12e", $u0->numify), "\n";
