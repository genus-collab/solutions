#!/usr/bin/perl
# t_curve_deform.pl -- dispersionless test: along the leaf curve, is the "deformed curve" u_t(lambda) of the closed form
#     u_t(lambda) = t (e^{t(lambda-1)} - 1)^tau (e^{t lambda} - 1)^{-1-tau} e^{-tau t (lambda + tau)}
# of the form  lambda + tau = Y_t(z) - sum_j c_j(t) z^j,  z = 1/u_t(lambda),  Y_t(z) = W(t) Y_0(z/W(t)),
# with W(t) = alpha(t) = -A_1(t)/c^2 the genus-0 scalar on the leaf curve, Y_0(z) = z - sigma sum_k Pi_k z^{-k}/(k+1)?
# Work at a numeric tau; expand in t to order M; at each t-order the coefficient is a Laurent series in 1/lambda (truncated).
# We use the variable mu = 1/lambda and represent z_t = 1/u_t as lambda * (series in t with coefficients polynomial in lambda) ...
# Simplest robust route: evaluate at many numeric lambda values (large) the truncated t-expansions and fit.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M) = @ARGV; $taus //= '2/7'; $M //= 3; my $tau = R($taus); my $sig = $tau*(1+$tau);
# --- leaf-curve data: P_j(t) and W(t)
my $x = -$tau/(1+$tau); my $c = 1 + $x;   # framing tau = -x/c
my $N = $M + 2;
my $Tt = T_of_t_series($N+1);
my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $N+1), $Tt, $N+1);
my %Pj;
for my $j (2..$N+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $N+1), $Tt, $N+1);
  my $p = [ R(1), map { R(0) } 1..$N+1 ]; $p = series_mul($p, $h1, $N+1) for 1..$j; my $rho = series_mul($hj, series_inv($p, $N+1), $N+1);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$N ]; }
# W solves W = 1 + sum_j Pi_j P_j W^j  (series in t)
my @W = (R(1), map { R(0) } 1..$N);
for my $it (1..$N+1) { my @nw = (R(1), map { R(0) } 1..$N); my $wp = [@W];
  for my $j (2..$N+1) { $wp = series_mul($wp, \@W, $N); my $wj = $wp; my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wj, $N); $nw[$_] = $nw[$_] + $q->[$_] for 0..$N; }
  @W = @nw; }
print "W(t) = ", join(', ', map { "$_" } @W[0..$M]), "\n";
# compare with alpha(t) = zeta(t) H_1(T)/t  (zeta = T e^{dphi/c})
my @u = (R(0), 4*$x/($c*$c), map { R(0) } 2..$N+1); my @omu = map { -$_ } @u; $omu[0] = R(1);
my $sq = series_pow(\@omu, R('1/2'), $N+1); my @den = map { $c * $_ } @$sq; $den[0] = $den[0] + $c;
my $ss = series_inv(\@den, $N+1); $ss = [ map { 2*$_ } @$ss ];
my @Ts = (R(0), @$ss[0..$N]); my @omTs = map { -$_ } @Ts; $omTs[0] = R(1);
my $dd = series_mul($ss, series_inv(\@omTs, $N+1), $N+1); my @dphi = (R(0)); $dphi[$_] = $dd->[$_-1]/$_ for 1..$N+1;
my $e = series_exp([ map { $_ / $c } @dphi ], $N+1); my @zT = (R(0), @$e[0..$N]); my $zeta = series_compose(\@zT, $Tt, $N+1);
my $al = series_mul([ map { $zeta->[$_+1] } 0..$N ], $h1, $N);
print "alpha(t) = ", join(', ', map { "$_" } @$al[0..$M]), "\n";
