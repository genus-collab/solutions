#!/usr/bin/perl
# r2_unstable.pl TAU D : the unstable cases of P2 (report s.5), independently:
#  (0,1): [x^mu] Int Y0 dX = sigma Pi_mu/(mu(mu+1))       (Y0 = -1/w, D = x d/dx, so [x^mu] F01 = [x^mu](-1/w)/mu)
#  (0,2): [x1^m1 x2^m2] log((w1-w2)/(x1-x2)) = sigma Pi_m1 Pi_m2/(m1+m2)   for m1, m2 >= 1, m1 + m2 <= D
# with w(x) from x = w (1+w)^(-tau-1) by Lagrange inversion ([x^m] w^k = (k/m) [w^(m-k)] (1+w)^(m(tau+1))).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my ($tau, $D) = @ARGV; $tau = R($tau); my $sig = $tau * (1 + $tau);
sub binR { my ($b, $i) = @_; my $r = R(1); for my $t (0..$i-1) { $r *= ($b - $t) / ($t + 1); } return $r; }
sub Pi_mu { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tau + $t) / $t; } return $r; }
my $N = $D + 2;
my @c = (R(0)); for my $m (1..$N) { push @c, binR($m * (1 + $tau), $m - 1) / $m; }        # w = Sum c_m x^m
# (0,1): -1/w = -(1/x) (1/(1 + c2 x + c3 x^2 + ...))
my @q = map { $c[$_ + 1] } 0..$N-1; my @inv = (R(1) / $q[0]);
for my $n (1..$N-1) { my $s = R(0); $s += $q[$_] * $inv[$n - $_] for 1..$n; push @inv, -$s / $q[0]; }
my ($b1, $c1) = (0, 0);
for my $mu (1..$D) { my $coef = -$inv[$mu + 1];   # [x^mu] (-1/w) = -[x^(mu+1)] (x/w)
  my $F = $coef / $mu; $c1++; $b1++ if $F != $sig * Pi_mu($mu) / ($mu * ($mu + 1)); }
# (0,2): S(x1,x2) = Sum_m c_m h_{m-1}(x1,x2) = 1 + ..., then log S
my %S; for my $m (1..$N) { for my $i (0..$m-1) { my $k = "$i," . ($m - 1 - $i); $S{$k} = ($S{$k} // R(0)) + $c[$m]; } }
sub bmul { my ($A, $B) = @_; my %C; for my $a (keys %$A) { my ($i1, $j1) = split /,/, $a; for my $b (keys %$B) { my ($i2, $j2) = split /,/, $b;
      my ($i, $j) = ($i1 + $i2, $j1 + $j2); next if $i + $j > $D; my $k = "$i,$j"; $C{$k} = ($C{$k} // R(0)) + $A->{$a} * $B->{$b}; } } return \%C; }
my %U = %S; delete $U{"0,0"}; for (keys %U) { delete $U{$_} if (split /,/)[0] + (split /,/)[1] > $D; }   # S = 1 + U
my %L; my $P = { "0,0" => R(1) };
for my $k (1..$D) { $P = bmul($P, \%U); for (keys %$P) { $L{$_} = ($L{$_} // R(0)) + $P->{$_} * R((-1) ** ($k + 1)) / $k; } }
my ($b2, $c2) = (0, 0);
for my $m1 (1..$D) { for my $m2 (1..$D - $m1) { my $v = $L{"$m1,$m2"} // R(0); $c2++; $b2++ if $v != $sig * Pi_mu($m1) * Pi_mu($m2) / ($m1 + $m2); } }
print "tau = $tau (sigma = $sig): (0,1) $c1 coefficients, $b1 mismatches; (0,2) $c2 coefficients (m1+m2 <= $D), $b2 mismatches\n";
