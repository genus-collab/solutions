#!/usr/bin/perl
# t_kp.pl -- the Marino-Vafa generating function R = sum_nu e^{tau x kappa_nu/2} s_nu(q^rho) s_nu(p) is a hypergeometric (content-product)
# tau function, hence a KP tau function in t_k = p_k/k.  Check of the first KP equation for C = log R at p = 0 and at the p_1 coefficient:
#   d1^4 C + 6 (d1^2 C)^2 + 12 d2^2 C - 12 d1 d3 C = 0     (d_k = d/dp_k).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my $tau = R('2/7'); my $P = 14;
# C_mu as x^{-|mu|} * series (mv_conn returns [p_mu] log R with valuation -|mu|)
sub Cm { my @mu = @_; return mv_conn(\@mu, $tau, $P); }
# constant term: 24 C_{1111} + 6 (2 C_{11})^2 + 12*2 C_{22} - 12 C_{31}; all have valuation -4 (products: -2-2)
my $t1 = ser_scal(R(24), Cm(1,1,1,1)); my $c11 = Cm(1,1); my $t2 = ser_scal(R(24), ser_mul($c11, $c11, $P));
my $t3 = ser_scal(R(24), Cm(2,2)); my $t4 = ser_scal(R(-12), Cm(3,1));
my $tot = ser_add(ser_add($t1, $t2), ser_add($t3, $t4));
print "KP at p=0 (coefficients of x^-4..x^", $P-4, "): ", join(' ', map { $tot->[$_] } 0..$P), "\n";
# coefficient of p_1: d/dp_1 of the KP expression at 0: 120 C_{1^5} + 6*2*(2C_{11})(6C_{111})... compute: d1 [ (d1^2 C)^2 ] = 2 (d1^2C)(d1^3C)
#   at p=0: d1^5 C = 120 C_{11111}; 2*(2 C_11)*(6 C_111) = 24 C_11 C_111; 12 d1 d2^2 C = 12*2 C_{221}; -12 d1^2 d3 C = -12*2 C_{311}
# valuations: C_{1^5}: -5; C_11 C_111: -2-3 = -5; C_221: -5; C_311: -5.
my $u1 = ser_scal(R(120), Cm(1,1,1,1,1)); my $u2 = ser_scal(R(6*24), ser_mul($c11, Cm(1,1,1), $P));
my $u3 = ser_scal(R(24), Cm(2,2,1)); my $u4 = ser_scal(R(-24), Cm(3,1,1));
my $tot2 = ser_add(ser_add($u1, $u2), ser_add($u3, $u4));
print "KP, p_1-coefficient: ", join(' ', map { $tot2->[$_] } 0..$P), "\n";
