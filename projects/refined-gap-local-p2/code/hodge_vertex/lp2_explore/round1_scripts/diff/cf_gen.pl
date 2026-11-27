#!/usr/bin/perl
# Exact conifold generators from the A-equation theta^2 A + 3y(3theta+1)(3theta+2)A = 0 in Del = 1+27y (theta = (Del-1) d/dDel):
# (N+1)^2 a_{N+1} = (2N^2+N+2/9) a_N - (N(N-1)+2/9) a_{N-1}, wc = sum a_n Del^n, a_0 = 1;  A_cf = -wc (BS v2: A_cf = -1 + O(xi)),
# t_cf' = wc/(1-Del) (so t_cf = Del + ...), F_cf = theta log A_cf, D_cf = X/A_cf^2, E_cf = (D_cf - X)/2, X = 1/Del.
use strict; use warnings; use Math::BigRat; my $N = 6; sub R { Math::BigRat->new(@_) }
my @a = (R(1)); my $am1 = R(0);
for my $n (0..$N) { my $prev = $n ? $a[$n-1] : R(0); $a[$n+1] = ((R(2*$n*$n+$n)+R(2,9))*$a[$n] - (R($n*($n-1))+R(2,9))*$prev)/R(($n+1)**2); }
print "wc coefficients a_0..a_4: ", join(", ", @a[0..4]), "\n";
# 1/(1-Del) * wc
my @b; for my $n (0..$N) { $b[$n] = R(0); $b[$n] += $a[$_] for 0..$n; }
my @t = (R(0)); $t[$_+1] = $b[$_]/R($_+1) for 0..$N-1;
print "t_cf coefficients of Del^1..Del^4: ", join(", ", @t[1..4]), "\n";
# F_cf(0) = theta log wc at 0 = (Del-1) wc'/wc |_0 = -a_1
print "F_cf(0) = XX_cf(0) = ", -$a[1], "\n";
# A_cf^2 = wc^2 ; 1/wc^2 series to Del^2
my @w2; for my $n (0..3) { $w2[$n]=R(0); $w2[$n] += $a[$_]*$a[$n-$_] for 0..$n; }
my @iw2 = (R(1)); for my $n (1..3) { my $s=R(0); $s += $w2[$_]*$iw2[$n-$_] for 1..$n; $iw2[$n] = -$s; }
print "D_cf = X/A_cf^2 = 1/Del + ", $iw2[1], " + ", $iw2[2], " Del + ...\n";
print "E_cf = (D_cf - X)/2 = ", $iw2[1]/R(2), " + ", $iw2[2]/R(2), " Del + ...\n";
