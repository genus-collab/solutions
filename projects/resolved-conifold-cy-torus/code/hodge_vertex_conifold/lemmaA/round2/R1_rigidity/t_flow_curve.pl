#!/usr/bin/perl
# t_flow_curve.pl -- test of the flow/curve conjecture (FC1): the flow phi_t of the lambda-line found from the ODE state satisfies
#     u_t(phi_t(lambda)) = u(lambda)/alpha(t),
#   u_t(lambda) = t (e^{t(lambda-1)} - 1)^tau (e^{t lambda} - 1)^{-1-tau} e^{-tau t (lambda+tau)}   (the curve of the closed form T1),
#   u(lambda) = (lambda-1)^tau lambda^{-1-tau},  alpha(t) = genus-0 scalar on the leaf curve.
# Equivalently  log u(phi) - log u(lambda) = K_t(phi) - log alpha(t),  K_t(mu) = (1+tau) Kappa(t mu) - tau Kappa(t(mu-1)).
# Solve for phi_t = lambda + sum_m t^m xi_m(lambda) order by order (exact, BigRat) and print xi_m.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M) = @ARGV; $taus //= '1'; $M //= 5; my $tau = R($taus); my $sig = $tau*(1+$tau);
# alpha(t) from the leaf curve (as in ode_state.pl)
my $Tt = T_of_t_series($M); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $M), $Tt, $M);
my %Pj; for my $j (2..$M+2) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $M), $Tt, $M); my $p = [R(1), map { R(0) } 1..$M]; $p = series_mul($p, $h1, $M) for 1..$j; my $rho = series_mul($hj, series_inv($p, $M), $M);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$M ]; }
my @al = (R(1), map { R(0) } 1..$M);
for my $it (0..$M+1) { my @nw = (R(1), map { R(0) } 1..$M); my $wp = [@al]; for my $j (2..$M+2) { $wp = series_mul($wp, \@al, $M); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp, $M); $nw[$_] = $nw[$_] + $q->[$_] for 0..$M; } @al = @nw; }
my $logal = series_log(\@al, $M);
# Bernoulli / kappa
my @Bn = (R(1)); for my $mm (1..$M+2) { my $s = R(0); for my $k (0..$mm-1) { $s = $s + binom($mm+1, $k) * $Bn[$k]; } $Bn[$mm] = -$s/($mm+1); }
my @kap = (R(0), $tau + R('1/2'), map { $Bn[$_]/($_*fact($_)) } 2..$M+2);
# polynomials in lambda as arrays (index = power).  phi = lambda + sum t^m xi_m.
sub padd { my ($a, $b, $s) = @_; $s //= R(1); my $n = @$a > @$b ? @$a : @$b; return [ map { ($a->[$_] // R(0)) + $s*($b->[$_] // R(0)) } 0..$n-1 ]; }
sub pmul { my ($a, $b) = @_; my @c = map { R(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i]->is_zero; for my $j (0..$#$b) { $c[$i+$j] = $c[$i+$j] + $a->[$i]*$b->[$j]; } } return \@c; }
sub pscal { my ($a, $s) = @_; return [ map { $_ * $s } @$a ]; }
# work with rational functions of lambda: equation in the form  sum_m t^m [ ... ] = 0 ; use the derivative form:
#   d/dt at each order.  Simplest: brute force with power series in t whose coefficients are rational functions represented as
#   (numerator polynomial)/(lambda(lambda-1))^e.  log u(phi) - log u(lambda) = tau log((phi-1)/(lambda-1)) - (1+tau) log(phi/lambda).
# With phi = lambda + d, d = O(t):  log(phi/lambda) = sum_k (-1)^{k+1} (d/lambda)^k/k,  log((phi-1)/(lambda-1)) = sum_k (-1)^{k+1} (d/(lambda-1))^k/k.
# Multiply through by (lambda(lambda-1))^M to stay polynomial: represent every t-coefficient as polynomial P/(lambda^a (lambda-1)^b) -- here we
# simply evaluate at several numeric lambda and fit; to keep it exact, solve at each order for xi_m via the linear term:
#   [t^m]: (u'/u)(lambda) xi_m + R_m(lambda) = [t^m](K_t(phi) - log alpha),  u'/u = -(lambda-1-tau)/(lambda(lambda-1)).
# All quantities are polynomials in lambda divided by powers of lambda and (lambda-1): use exact rational-function arithmetic via evaluation
# at sample points lambda_s and polynomial interpolation of xi_m (degree <= m+1).
my @lams = map { R($_) + R(1)/R(7) } 2..($M+8);
my %xi;   # xi{m}[s] values at sample points
for my $m (1..$M) {
  my @vals;
  for my $s (0..$#lams) { my $L = $lams[$s];
    # d(t) = sum_{k<m} t^k xi_k(L) (known), as series in t
    my @d = (R(0)); for my $k (1..$M) { $d[$k] = ($k < $m) ? $xi{$k}[$s] : R(0); }
    # LHS series: tau log(1 + d/(L-1)) - (1+tau) log(1 + d/L)
    my @a1 = map { $_ / ($L - 1) } @d; $a1[0] = R(1); my @a2 = map { $_ / $L } @d; $a2[0] = R(1);
    my $l1 = series_log(\@a1, $M); my $l2 = series_log(\@a2, $M);
    my @lhs = map { $tau*$l1->[$_] - (1+$tau)*$l2->[$_] } 0..$M;
    # RHS: K_t(phi) - log alpha, with phi = L + d:  K_t(phi) = sum_k kappa_k t^k [ (1+tau) phi^k - tau (phi-1)^k ]
    my @phi = ($L, @d[1..$M]); my @phim1 = ($L - 1, @d[1..$M]);
    my @rhs = map { -$logal->[$_] } 0..$M; $rhs[0] = R(0);
    my $pp = [R(1), map { R(0) } 1..$M]; my $pq = [R(1), map { R(0) } 1..$M];
    for my $k (1..$M) { $pp = series_mul($pp, \@phi, $M); $pq = series_mul($pq, \@phim1, $M);
      for my $i (0..$M-$k) { $rhs[$i+$k] = $rhs[$i+$k] + $kap[$k]*((1+$tau)*$pp->[$i] - $tau*$pq->[$i]); } }
    # linear coefficient of xi_m in lhs at t^m: (u'/u)(L) = -(L-1-tau)/(L(L-1))
    my $up = -($L - 1 - $tau)/($L*($L - 1));
    push @vals, ($rhs[$m] - $lhs[$m]) / $up; }
  $xi{$m} = \@vals;
  # interpolate as polynomial of degree <= m+1 using the first m+2 points, check on the rest
  my $deg = $m + 1; my @c = interp([ @lams[0..$deg] ], [ @vals[0..$deg] ]);
  my $ok = 1; for my $s ($deg+1..$#lams) { my $v = R(0); my $pw = R(1); for my $cc (@c) { $v = $v + $cc*$pw; $pw = $pw*$lams[$s]; } $ok = 0 if $v != $vals[$s]; }
  print "xi_$m = ", join(' + ', map { "($c[$_])l^$_" } 0..$deg), "   (polynomial of degree <= ", $m+1, ": ", ($ok ? "yes" : "NO"), ")\n"; }
sub interp { my ($x, $y) = @_; my $n = @$x; my @c = map { R(0) } 0..$n-1;
  for my $i (0..$n-1) { my @num = (R(1)); my $den = R(1); for my $j (0..$n-1) { next if $j == $i; @num = @{ pmul(\@num, [ -$x->[$j], R(1) ]) }; $den = $den * ($x->[$i] - $x->[$j]); }
    for my $k (0..$#num) { $c[$k] = $c[$k] + $y->[$i] * $num[$k] / $den; } } return @c; }
