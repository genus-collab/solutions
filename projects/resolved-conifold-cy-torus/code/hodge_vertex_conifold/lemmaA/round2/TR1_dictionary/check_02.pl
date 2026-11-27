#!/usr/bin/perl
# check_02.pl tau [D]: the unstable (0,2) case of the dictionary.
#   F_{0,2}(x1,x2) := log( (w(x1) - w(x2)) / (x1 - x2) ) - log w'(0)   (so d1 d2 F_{0,2} = B - dx1 dx2/(x1-x2)^2),
#   w(x) the inverse of x = w/(1+w)^(tau+1) at w = 0.  Claim: [x1^m1 x2^m2] F_{0,2} = sigma Pi_m1 Pi_m2 / (m1 + m2)
#   = sigma^(1-g) Pi Pi Xi_0(m1,m2) with g = 0, Xi_0(m1,m2) = 1/(m1+m2) (Route A's Gamma(J)/Gamma(J-r+3) at r = 2).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my ($tau, $D) = @ARGV; $tau = R($tau // 1); $D //= 8;
my $sigma = $tau * ($tau + 1);
sub binom_r { my ($N, $i) = @_; my $r = R(1); for my $t (1..$i) { $r *= ($N - $t + 1) / $t; } return $r; }
# Lagrange: [x^m] w = (1/m) [w^(m-1)] (1+w)^(m(tau+1)) = binom(m(tau+1), m-1)/m
my @c = (R(0)); for my $m (1..$D+1) { push @c, binom_r($m * ($tau + 1), $m - 1) / $m; }
# bivariate series P(x1,x2) = Sum_m c_m h_{m-1}(x1,x2), total degree <= D ; key "i,j"
my %P; for my $m (1..$D+1) { for my $i (0..$m-1) { my $j = $m - 1 - $i; $P{"$i,$j"} = ($P{"$i,$j"} // R(0)) + $c[$m]; } }
# log P with P(0,0) = 1: log(1+Q) = Sum (-1)^(k-1) Q^k / k
my %Q = %P; $Q{"0,0"} = R(0);
sub bmul { my ($A, $B) = @_; my %C; for my $ka (keys %$A) { my ($i1,$j1) = split /,/, $ka; for my $kb (keys %$B) { my ($i2,$j2) = split /,/, $kb;
      next if $i1 + $i2 + $j1 + $j2 > $D; my $k = ($i1+$i2) . "," . ($j1+$j2); $C{$k} = ($C{$k} // R(0)) + $A->{$ka} * $B->{$kb}; } } return \%C; }
my %F; my $pw = { "0,0" => R(1) };
for my $k (1..$D) { $pw = bmul($pw, \%Q); for my $key (keys %$pw) { $F{$key} = ($F{$key} // R(0)) + R((-1) ** ($k - 1)) / $k * $pw->{$key}; } }
my ($n, $bad) = (0, 0);
sub Pi { my $mu = shift; my $r = R(1); for my $t (1..$mu-1) { $r *= ($mu * $tau + $t) / $t; } return $r; }
for my $m1 (1..$D-1) { for my $m2 (1..$D - $m1) { my $got = $F{"$m1,$m2"} // R(0); my $exp = $sigma * Pi($m1) * Pi($m2) / ($m1 + $m2);
    $n++; $bad++ if $got != $exp; } }
print "tau = $tau: (0,2) dictionary  [x1^m1 x2^m2] F_{0,2} = sigma Pi_m1 Pi_m2/(m1+m2): $n coefficients (m1+m2 <= $D), mismatches $bad\n";
