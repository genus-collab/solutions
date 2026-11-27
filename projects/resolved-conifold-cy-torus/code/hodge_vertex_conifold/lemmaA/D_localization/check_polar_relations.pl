# check_polar_relations.pl -- algebraic checks behind report.md section 5 (polar relations T_k):
#  (1) the explicit leaf coefficients a_j^(k) of rho q_t (leaf_structure.md, k <= 3; leaf_structure_proof.md Cor. 11.3) satisfy
#      top:      a_{k+1}^(k) = (-1)^(k+1) C_{k+1},  C_j = (c^2/j) eps_j(sigma),  eps_j = C(j tau + j - 1, j - 1)^2 = prod_{k<j} (j^2 sigma + k(j-k)) / ((j-1)!)^2
#      sub-top:  a_k^(k) = a_k^(k-1) * k (1-x)^2/(2c^2)  and  a_1^(1) = -(1-x)^2/2   (first-order string flow s = (1-x)^2 t/(2c))
#  (2) in [t^k] F_h(rho q_t) (c-normalisation, string equation Xi(1,R) = (J_R/c^2) Xi(R)) the coefficients of Xi(R) with D(R) = k-1 vanish
#      identically for k = 2, 3 (Xi(2) at k = 2; Xi(3), Xi(2,2) at k = 3), and the residual Xi(2)-coefficient at k = 3 is -(1-x)^2 x^2/(12 c^4);
#  (3) the c = 1 forms of T_2 and T_3 used in the report.
# Evaluated exactly at several rational x (tau = -x/c, sigma = tau(1+tau)).
use strict; use warnings; use Math::BigRat;
sub Q { Math::BigRat->new(@_) }
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; return $r; }
my $ok = 1;
sub chk { my ($name, $a, $b) = @_; my $good = ($a == $b); $ok &&= $good; printf "%-60s %s\n", $name, $good ? "ok" : "MISMATCH ($a vs $b)"; }
for my $x (Q(2), Q(1,3), Q(-5,7), Q(3,2)) {
  my $c = 1 + $x; my $tau = -$x / $c; my $sig = $tau * (1 + $tau);
  print "x = $x (c = $c, tau = $tau, sigma = $sig)\n";
  # explicit coefficients from leaf_structure.md (exact fits verified there), c-normalisation
  my %a;
  $a{1}{1} = -(1-$x)**2/2;             $a{2}{1} = (1-$x)**2/2;
  $a{1}{2} = -(2 - 4*$x + 3*$x**2 - 4*$x**3 + 2*$x**4)/(12*$c**2);
  $a{2}{2} = 6*(1-$x)**4/(12*$c**2);   $a{3}{2} = -(2-$x)**2*(1-2*$x)**2/(12*$c**2);
  my $p3 = (1-$x)**2/$c**4;
  $a{1}{3} = -$p3*(3 + 5*$x**2 + 3*$x**4)/72;
  $a{2}{3} = $p3*(7 - 20*$x + 22*$x**2 - 20*$x**3 + 7*$x**4)/24;
  $a{3}{3} = -$p3*(2-$x)**2*(1-2*$x)**2/8;
  $a{4}{3} = $p3*(3-$x)**2*(1-3*$x)**2/36;
  # sum rule sum_j a_j^(k) = 0 for k >= 1 (Cor. 11.2(a))
  for my $k (1..3) { my $s = 0; $s += $a{$_}{$k} // 0 for 1..4; chk("  sum rule k=$k", $s, Q(0)); }
  # top coefficients
  my $eps = sub { my $j = shift; my $p = Q(1); for my $k (1..$j-1) { $p *= ($j*$j*$sig + $k*($j-$k)); } return $p / fact($j-1)**2; };
  my $C = sub { my $j = shift; return $c**2/$j * $eps->($j); };
  for my $k (1..3) { chk("  top a_{".($k+1)."}^($k) = (-1)^".($k+1)." C_".($k+1), $a{$k+1}{$k}, Q((-1)**($k+1)) * $C->($k+1)); }
  # sub-top
  chk("  a_1^(1) = -(1-x)^2/2 = -c s1", $a{1}{1}, -(1-$x)**2/2);
  for my $k (2..3) { chk("  sub-top a_$k^($k) = a_$k^(".($k-1).") k (1-x)^2/(2c^2)", $a{$k}{$k}, $a{$k}{$k-1} * $k * (1-$x)**2/(2*$c**2)); }
  # (2) coefficients in [t^k] F(rho q_t); Xi(1,R) = (J/c^2) Xi(R)
  my $J2 = 2/$c**2;  # Xi(1,2) = (2/c^2) Xi(2)
  chk("  k=2: coefficient of Xi(2) vanishes", $a{2}{2} + $a{1}{1}*$a{2}{1}*$J2, Q(0));
  chk("  k=3: coefficient of Xi(3) vanishes", $a{3}{3} + 3/$c**2 * $a{1}{1}*$a{3}{2}, Q(0));
  chk("  k=3: coefficient of Xi(2,2) vanishes", $a{2}{1}*$a{2}{2} + 2*$a{1}{1}*$a{2}{1}**2 * 4/(4*$c**2) * 1, Q(0));
  # (Xi(1,2,2) = (4/c^2) Xi(2,2); term (1/6) * 3 a1 a2^2 * (4/c^2) = 2 a1 a2^2/c^2 )
  my $res2 = $a{2}{3} + 2*($a{1}{1}*$a{2}{2} + $a{2}{1}*$a{1}{2})/$c**2 + 3*$a{1}{1}**2*$a{2}{1}/$c**4;
  chk("  k=3: residual Xi(2) coefficient = -(1-x)^2 x^2/(12 c^4)", $res2, -(1-$x)**2*$x**2/(12*$c**4));
  # (3) c = 1 forms: T_2: -C_3 Xi^c(3) + C_2^2/2 Xi^c(2,2) = 0  <=>  2(2+9s)^2 Xi(3) = 3(1+4s)^2 Xi(2,2)  (Xi^c(R) = c^-(3h-3+2r) Xi(R))
  #     check proportionality of the coefficient vectors after the c-powers: (-C_3 c^-2, C_2^2/2 c^-4) ~ (-2(2+9s)^2, 3(1+4s)^2)
  my @v = (-$C->(3)/$c**2, $C->(2)**2/2/$c**4); my @w = (-2*(2+9*$sig)**2, 3*(1+4*$sig)**2);
  chk("  T_2 c=1 form proportional", $v[0]*$w[1], $v[1]*$w[0]);
  # T_3: C_4 Xi^c(4) - C_2 C_3 Xi^c(2,3) + C_2^3/6 Xi^c(2,2,2): c-powers c^-2, c^-4, c^-6  ~ (4(3+16s)^2, -6(2+9s)^2, 3(1+4s)^2) times (1+4s)/144
  my @v3 = ($C->(4)/$c**2, -$C->(2)*$C->(3)/$c**4, $C->(2)**3/6/$c**6);
  my @w3 = (4*(3+16*$sig)**2, -6*(2+9*$sig)**2, 3*(1+4*$sig)**2);
  chk("  T_3 c=1 form proportional (1,2)", $v3[0]*$w3[1], $v3[1]*$w3[0]);
  chk("  T_3 c=1 form proportional (1,3)", $v3[0]*$w3[2], $v3[2]*$w3[0]);
}
print $ok ? "ALL OK\n" : "SOME MISMATCH\n";
