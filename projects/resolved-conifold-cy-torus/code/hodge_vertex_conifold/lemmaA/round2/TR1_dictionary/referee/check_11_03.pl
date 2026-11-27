#!/usr/bin/perl
# Referee check of the base cases (0,3), (1,1) at many tau (independent code TRw.pm), against the closed forms
#   Y_0 = -1/w      : omega_{0,3} = sigma dxi_0^3,   omega_{1,1} = -((2+1/sigma)/24) dxi_0 + (1/24) dxi_1
#   y_B = log(1+w)  : omega_{0,3} = sigma^2 dxi_0^3, omega_{1,1} = -((1+sigma)/24) dxi_0 + (sigma/24) dxi_1
# and against the general one-branch-point formula derived by the referee:
#   omega_{1,1} = (1/(4 c2 y1)) [ (beta0 - y3/(4 y1)) e_0 + e_1/12 ],  omega_{0,3} = e_0^3/(2 c2 y1),
#   (zeta-coordinate data: X - X(a) = -c2 zeta^2, y = y(a) + y1 zeta + y2 zeta^2 + y3 zeta^3 + ..., beta0 = u3 - u2^2).
# Comparison through x-coefficients: [x^mu] F_{1,1} = C0 Pi_mu + C1 mu Pi_mu (mu = 1, 2, 3, 4), [x^mu] F_{0,3} = C Pi Pi Pi.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use TRw;
sub R { Math::BigRat->new(@_) }
my @taus = qw(1 2 3 4 7 1/2 1/3 2/3 3/5 5/2 7/3 13/5 -1/2 -1/3 -2/5 -2 -3 -5/2 -7/3 -1/7 -6/7 11 -11/3 9/4);
my ($n, $bad) = (0, 0);
for my $t (@taus) { for my $y ('base', 'bkmp') {
    TRw::setup(tau => $t, L => 24, y => $y); TRw::initial(); my $sig = $TRw::sigma;
    my ($C0, $C1, $C03) = $y eq 'base' ? (-(2 + 1/$sig)/24, R(1)/24, $sig) : (-(1 + $sig)/24, $sig/24, $sig ** 2);
    for my $mu (1..4) { my $pi = TRw::Pi_mu($mu); my $tr = TRw::xcoef($TRw::W{"1,1"}, $mu); my $ex = ($C0 + $C1 * $mu) * $pi; $n++; $bad++ if $tr != $ex; }
    for my $m ([1,1,1], [1,2,1], [2,2,1], [3,1,2]) { my $pi = R(1); $pi *= TRw::Pi_mu($_) for @$m; my $tr = TRw::xcoef($TRw::W{"0,3"}, @$m); $n++; $bad++ if $tr != $C03 * $pi; } } }
print "base cases (0,3), (1,1): ", scalar(@taus), " values of tau x 2 curves, $n x-coefficients compared, mismatches $bad\n";
