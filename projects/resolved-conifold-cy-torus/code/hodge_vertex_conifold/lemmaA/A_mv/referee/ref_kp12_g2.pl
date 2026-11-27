#!/usr/bin/perl
# ref_kp12_g2.pl -- referee check (independent of the route's code): the first two KP equations, evaluated at the dilaton point
# at genus G = 2 with R = empty, using ONLY the exact Hodge tables lemmaA/data/q_g2_r1.txt, q_g2_r2.txt and genus-0/1 values
# derived by hand in the referee report:
#   Xi_0(1,1) = 1/2,  Xi_1(1) = -(1+1/sigma)/24,  Xi_1(2) = -1/(24 sigma)     (Q_1(s;(j)) = (j - 1 - 1/sigma + s)/24)
# KP1 (p-derivatives): F_1111 + 6 F_11^2 + 12 F_22 - 12 F_13 = 0  ->
#   X1 := sigma [6 Xi_1(1) + 6 Xi_1(1)^2] + 12 Pi_2^2 Xi_2(2,2) - 36 Pi_3 Xi_2(3)                         = 0 ?
# KP2 (p-derivatives): F_1112 + 6 F_11 F_12 + 6 F_23 - 6 F_14 = 0  ->
#   X2 := sigma Pi_2 [24 Xi_1(2) + 6 (Xi_2(2) + 2 Xi_1(1) Xi_1(2))] + 6 Pi_2 Pi_3 Xi_2(2,3) - 24 Pi_4 Xi_2(4)  = 0 ?
# c-normalisation of the tables: Xi^(c=1) = c^(3h-3+2r) Xi^(c-norm), c = 1+x = 1/(1+tau), x = -tau/(1+tau).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rparse; use Math::BigRat;
my $data = "$FindBin::Bin/../../data";
my ($q1) = read_table("$data/q_g2_r1.txt"); my ($q2) = read_table("$data/q_g2_r2.txt");
die "header" unless $q1->[0] == 2 && $q1->[1] == 1 && $q2->[0] == 2 && $q2->[1] == 2;
sub R { Math::BigRat->new(@_) }
for my $ts ('2/7', '3/2', '5', '-7/3', '1/9', '-4/11') {
  my $tau = R($ts); my $x = -$tau / (1 + $tau); my $c = 1 + $x; my $sig = $tau * (1 + $tau);
  my $Pi = sub { my $j = shift; my $v = R(1); for my $b (1..$j-1) { $v *= ($j * $tau + $b); } for my $b (1..$j-1) { $v /= $b; } return $v; };
  my $X1v = sub { my $j = shift; return $c**5 * eval_expr($q1->[2], { s => -1, j1 => $j, x => "$x" }); };
  my $X2v = sub { my ($a, $b) = @_; return $c**7 * eval_expr($q2->[2], { s => -1, j1 => $a, j2 => $b, x => "$x" }); };
  my $xi11 = -(1 + 1/$sig) / 24; my $xi12 = -1 / (24 * $sig);
  my ($x2, $x3, $x4, $x22, $x23) = ($X1v->(2), $X1v->(3), $X1v->(4), $X2v->(2,2), $X2v->(2,3));
  my $K1 = $sig * (6*$xi11 + 6*$xi11**2) + 12 * $Pi->(2)**2 * $x22 - 36 * $Pi->(3) * $x3;
  my $K2 = $sig * $Pi->(2) * (24*$xi12 + 6*($x2 + 2*$xi11*$xi12)) + 6 * $Pi->(2) * $Pi->(3) * $x23 - 24 * $Pi->(4) * $x4;
  # sensitivity: the same with a wrong genus-1 value
  my $K1w = $sig * (6*$xi11*R('1.5') + 6*$xi11**2) + 12 * $Pi->(2)**2 * $x22 - 36 * $Pi->(3) * $x3;
  printf "tau=%-6s Xi_2(2)=%s  Xi_2(3)=%s  Xi_2(2,2)=%s\n    KP1 total = %s   KP2 total = %s   (KP1 with Xi_1(1) scaled by 3/2: %s)\n",
    $ts, $x2, $x3, $x22, $K1, $K2, $K1w;
}
