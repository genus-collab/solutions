#!/usr/bin/perl
# Critic spot check: rank of the 2x2 block E_11(a, a'), a, a' in {0,1}, near t = 0.
# Input: the fitted leaf coefficients of [t^k](rho f_{1,b}) for b = 0, 1, k = 0, 1, 2
# (hvc_leaf_structure_data.m, Mathematica fits from the stored code factors).
# [psi^a] lam_j = j^a / c^(a+1).  det(rho E-block) = rho^2 det(E-block); since det = O(t^2),
# its t^2 coefficient equals that of det(E-block) and of det(rho^2 E-block).
# Prediction (from the loop table e^(k)): [t^2] det = (1-x)^4 x^6 / (18 c^10).
use strict; use warnings;
use Math::BigRat;
sub R { Math::BigRat->new(@_) }
sub pw { my ($b,$e)=@_; my $r=R(1); $r*=$b for 1..$e; return $r; }
for my $xs ('2','1/3','3/7','-5/2','1') {
  my $x = R($xs); my $c = $x + 1;
  my $om = 1 - $x;
  # leaf coefficients: $L{b}{k} = [coef lam_1, coef lam_2, coef lam_3]
  my %L;
  $L{0}{0} = [ -pw($x,3)/(2*$c) ];
  $L{0}{1} = [ -pw($om,2)*pw($x,3)/(12*pw($c,3)), 2*pw($om,2)*pw($x,3)/(3*pw($c,3)) ];
  $L{0}{2} = [ pw($x,4)*(4-3*$x+4*$x*$x)/(48*pw($c,5)), pw($om,4)*pw($x,3)/(3*pw($c,5)),
               -3*pw(2-$x,2)*pw($x,3)*pw(1-2*$x,2)/(16*pw($c,5)) ];
  $L{1}{0} = [ -pw($x,3)/(2*pw($c,2)) ];
  $L{1}{1} = [ 7*pw($om,2)*pw($x,3)/(12*pw($c,4)), 2*pw($om,2)*pw($x,3)/(3*pw($c,4)) ];
  $L{1}{2} = [ -pw($x,3)*(8-108*$x+213*$x*$x-108*pw($x,3)+8*pw($x,4))/(48*pw($c,6)),
               -2*pw($om,4)*pw($x,3)/(3*pw($c,6)), -3*pw(2-$x,2)*pw($x,3)*pw(1-2*$x,2)/(16*pw($c,6)) ];
  # M[k][a][b] = [psi^a][t^k](rho f_{1,b})
  my @M;
  for my $k (0..2) { for my $a (0..1) { for my $b (0..1) {
    my $s = R(0); my $co = $L{$b}{$k};
    for my $j (1..scalar(@$co)) { $s += $co->[$j-1] * pw(R($j),$a) / pw($c,$a+1); }
    $M[$k][$a][$b] = $s; } } }
  # symmetry check E(a,b) = E(b,a) at each order (rho f: same factor rho on both, so symmetric)
  my $sym = 1; for my $k (0..2) { $sym = 0 unless $M[$k][0][1] == $M[$k][1][0]; }
  # det series to t^2
  my @d = (R(0),R(0),R(0));
  for my $i (0..2) { for my $j (0..2) { next if $i+$j>2;
     $d[$i+$j] += $M[$i][0][0]*$M[$j][1][1] - $M[$i][0][1]*$M[$j][1][0]; } }
  my $pred = pw($om,4)*pw($x,6)/(18*pw($c,10));
  printf "x=%-5s sym=%d  det t^0=%s  t^1=%s  t^2=%s  pred=%s  %s\n", $xs, $sym, $d[0], $d[1], $d[2], $pred,
     ($d[2]==$pred ? "MATCH" : "MISMATCH");
}
