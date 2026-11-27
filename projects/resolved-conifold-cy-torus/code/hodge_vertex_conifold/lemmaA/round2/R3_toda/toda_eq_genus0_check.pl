#!/usr/bin/perl
# Line (2): genus-0, empty-profile part of the CY Toda equation d_{p1}(e^{-H} Eq e^H) = exp(H(e^b o p) - 2H + H(e^{-b} o p)):
#   sigma d/dy [w^0] exp(-Y(w)) = exp( (tau/(1+tau)) T(y) ),
#   Y(w) = w y/sigma + Sum_k k Pi_k w^{-k} G0_k(y),  G0_k(y) = Sum_a (k+a)^(a-2) y^a/a!,  T(y) = Sum_a a^(a-1) y^a/a! (tree function).
# Exact check as power series in y through y^M.  (At the dilaton point T = u, so the right side carries e^{tau u/(1+tau)}.)
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $tau = R($ARGV[0] // "2/7"); my $M = $ARGV[1] // 8;
my $sig = $tau * ($tau + 1);
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub Pi { my $k = shift; my $p = R(1); for my $a (1..$k-1) { $p *= ($k*$tau + $a); } $p / fact($k-1) }
sub smul { my ($a, $b) = @_; my @c = map { R(0) } 0..$M; for my $i (0..$M) { next if $a->[$i] == 0; for my $j (0..$M-$i) { $c[$i+$j] += $a->[$i] * $b->[$j]; } } \@c }
sub sexp { my ($e) = @_; my @f = (R(1)); for my $n (1..$M) { my $s = R(0); for my $m (1..$n) { $s += $m * $e->[$m] * $f[$n-$m]; } $f[$n] = $s / $n; } [ map { $_ // R(0) } @f ] }   # needs e[0] = 0
my $K = $M + 1;   # w^{-k}, k <= K suffices: [w^0] pairs w^m (from e^{-wy/sigma}, order y^m) with w^{-m}
my @G; for my $k (1..$K) { $G[$k] = [ map { my $a = $_; R($k + $a) ** ($a - 2) / fact($a) } 0..$M ]; }
# c_m(y) = [w^{-m}] exp(-Sum_k k Pi_k w^{-k} G_k(y)); compute by the recursion for exp of a series in v = 1/w:
#   E(v) = exp(L(v)), L(v) = -Sum_k k Pi_k G_k v^k;  m c_m = Sum_{k=1}^m k L_k c_{m-k}
my @L; $L[$_] = [ map { -$_ * $_ } () ] for (); for my $k (1..$K) { my $c = -$k * Pi($k); $L[$k] = [ map { $_ * $c } @{$G[$k]} ]; }
my @c; $c[0] = [ R(1), map { R(0) } 1..$M ];
for my $m (1..$M) { my $s = [ map { R(0) } 0..$M ]; for my $k (1..$m) { my $t = smul($L[$k], $c[$m-$k]); $s->[$_] += $k * $t->[$_] for 0..$M; } $c[$m] = [ map { $_ / $m } @$s ]; }
# [w^0] exp(-Y) = Sum_m (-y/sigma)^m/m! c_m(y)
my @W = map { R(0) } 0..$M;
for my $m (0..$M) { my $f = (R(-1)/$sig) ** $m / fact($m); for my $i (0..$M-$m) { $W[$i+$m] += $f * $c[$m][$i]; } }
my @lhs = map { R(0) } 0..$M-1; for my $i (1..$M) { $lhs[$i-1] = $sig * $i * $W[$i]; }
my @T = map { $_ == 0 ? R(0) : R($_) ** ($_ - 1) / fact($_) } 0..$M;
my $rhs = sexp([ map { $_ * $tau / (1 + $tau) } @T ]);
my $ok = 1; for my $i (0..$M-1) { $ok = 0 if $lhs[$i] != $rhs->[$i]; }
print "lhs: @lhs[0..4] ...\nrhs: @{$rhs}[0..4] ...\n";
printf "genus-0 empty-profile CY Toda equation at tau=%s through y^%d: %s\n", $tau, $M-1, ($ok ? "OK" : "FAIL");
