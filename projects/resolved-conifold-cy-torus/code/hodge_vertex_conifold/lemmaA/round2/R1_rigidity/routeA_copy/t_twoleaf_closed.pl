#!/usr/bin/perl
# t_twoleaf_closed.pl -- algebraic check of the closed form of the two-leaf KP relation (R = empty) used in report section 3.6:
#   sum_{a<=al,b<=be} Pi_a Pi_b Xi(a,b) E0(al-a,be-b) = (al+be-1) Pi_{al+be-1} Xi(al+be-1)       (Theorem D, R = empty)
# is equivalent to    A(z,w) := sum_{a,b>=1} Pi_a Pi_b Xi(a,b) z^-a w^-b = (V(w) - V(z)) / (Y(z) - Y(w)),
#   V(z) = sum_c c Pi_c Xi(c) z^-c,  Y(z) = z - U(z),  U(z) = sigma sum_k Pi_k z^-k/(k+1),   Xi(1,b) = b Xi(b), Xi(1,1) = 0.
# Test with random rational one-leaf values: solve the relations recursively for Xi(a,b), compare with the series expansion of the closed form.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($tau, $M) = (R('3/5'), 8); my $sig = $tau*(1+$tau);
sub Pi { my $j = shift; return R(1) if $j <= 1; my $p = R(1); for my $b (1..$j-1) { $p *= ($j*$tau + $b); } return $p / fact($j-1); }
srand(7); my %X1; $X1{1} = R(0); for my $c (2..2*$M) { $X1{$c} = R(int(rand(19)) - 9) / R(1 + int(rand(7))); }
# E0(m,n) = [z^-m w^-n] exp(sum sigma Pi_a Pi_b z^-a w^-b/(a+b))
my %Om; for my $a (1..$M) { for my $b (1..$M) { $Om{"$a,$b"} = $sig * Pi($a) * Pi($b) / ($a + $b); } }
my %E = ("0,0" => R(1)); my %P = ("0,0" => R(1));
for my $n (1..$M) { my %Q; for my $k1 (keys %P) { my ($a1,$b1) = split /,/, $k1; for my $k2 (keys %Om) { my ($a2,$b2) = split /,/, $k2; my ($a,$b) = ($a1+$a2,$b1+$b2); next if $a > $M || $b > $M; $Q{"$a,$b"} = ($Q{"$a,$b"} // R(0)) + $P{$k1}*$Om{$k2}; } } %P = %Q; for my $k (keys %P) { $E{$k} = ($E{$k} // R(0)) + $P{$k}/fact($n); } }
sub E0 { my ($m,$n) = @_; return $E{"$m,$n"} // R(0); }
# solve recursively for X2(a,b) (a,b >= 2), with X2(1,b) = b X1(b)
my %X2; sub X2 { my ($a,$b) = @_; return R(0) if $a == 1 && $b == 1; return $b * $X1{$b} if $a == 1; return $a * $X1{$a} if $b == 1; return $X2{"$a,$b"}; }
for my $s (4..2*$M) { for my $a (2..$s-2) { my $b = $s - $a; next if $a > $M || $b > $M;
    my $rhs = ($a+$b-1) * Pi($a+$b-1) * $X1{$a+$b-1}; my $lhs = R(0);
    for my $a1 (1..$a) { for my $b1 (1..$b) { next if $a1 == $a && $b1 == $b; my $v = X2($a1,$b1); next unless defined $v; $lhs += Pi($a1)*Pi($b1)*$v*E0($a-$a1,$b-$b1); } }
    $X2{"$a,$b"} = ($rhs - $lhs) / (Pi($a)*Pi($b)); } }
# closed form: (V(w)-V(z))/(Y(z)-Y(w)) in u = 1/z, v = 1/w : Y(z)-Y(w) = (z - w) - (U(z)-U(w)) ; V(w)-V(z) = sum_c cPi_c X1(c)(v^c - u^c)
# write z - w = (v - u)/(uv); divide numerator and denominator by (v - u): N = sum_c cPi_c X1(c) h_c(u,v), Dn = 1/(uv) - sum_k sigma Pi_k/(k+1) h_k(u,v)
# with h_c = sum_{i+j=c-1} u^i v^j:  B = uv sum_c c Pi_c X1(c) h_c,  E = 1 + uv sum_k sigma Pi_k/(k+1) h_k,  A = B/E.
my %Nn; for my $c (1..2*$M) { my $cc = $c * Pi($c) * $X1{$c}; for my $i (0..$c-1) { my $j = $c-1-$i; next if $i+1 > $M || $j+1 > $M; $Nn{($i+1).",".($j+1)} = ($Nn{($i+1).",".($j+1)} // R(0)) + $cc; } }
my %Dg; for my $k (1..2*$M) { my $cc = $sig * Pi($k) / ($k+1); for my $i (0..$k-1) { my $j = $k-1-$i; next if $i+1 > $M || $j+1 > $M; $Dg{($i+1).",".($j+1)} = ($Dg{($i+1).",".($j+1)} // R(0)) + $cc; } }
# A = Nn / E = Nn / (1 + Dg) = Nn * sum_n (-Dg)^n   (E = 1 + (U(z)-U(w))/(w-z) = 1 + Dg)
my %A = %Nn; my %T = %Nn;
for my $n (1..$M) { my %Q; for my $k1 (keys %T) { my ($a1,$b1) = split /,/, $k1; for my $k2 (keys %Dg) { my ($a2,$b2) = split /,/, $k2; my ($a,$b) = ($a1+$a2,$b1+$b2); next if $a > $M || $b > $M; $Q{"$a,$b"} = ($Q{"$a,$b"} // R(0)) + $T{$k1}*$Dg{$k2}; } } %T = %Q; for my $k (keys %T) { $A{$k} = ($A{$k} // R(0)) + R((-1)**$n) * $T{$k}; } }
my ($n, $bad) = (0, 0);
for my $a (1..$M) { for my $b (1..$M) { next if $a + $b > $M + 1; my $lhs = Pi($a)*Pi($b)*(X2($a,$b) // R(0)); my $rhs = $A{"$a,$b"} // R(0); $n++; $bad++ if $lhs != $rhs; } }
print "closed form of the two-leaf KP relation: $n coefficients compared, $bad mismatches\n";
