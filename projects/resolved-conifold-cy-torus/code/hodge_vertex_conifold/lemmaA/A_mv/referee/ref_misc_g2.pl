#!/usr/bin/perl
# ref_misc_g2.pl -- referee spot checks from the exact Hodge tables only (no route code):
#  (a) Theorem B (2-insertion relation) at genus 2, R = (3) and R = (2), outside the range (all terms non-zero):
#      R=(3): Xi_2(3,2) - 4 Pi_4/(Pi_3 s') Xi_2(4) + (2 sigma/Pi_3) Xi_1(2) (2 + Xi_1(1)) = 0
#      R=(2): Xi_2(2,2) - 3 Pi_3/s'^2 Xi_2(3) + sigma/(2 s'^2) Xi_1(1)(1 + Xi_1(1)) = 0      (s' = 2 tau + 1)
#  (b) simple zero of Q_h(s;R) at s = -1 for the irreducible points of genus 2, 3 (exact, tables q_g2_r*, q_g3_r*)
#  (c) the lambda_1-component Xi^(1)_2(j1,j2) at (1/2, 5/2) (claimed -1/3072 by the route)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rparse; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $data = "$FindBin::Bin/../../data"; my %T; for my $f (qw(q_g2_r1 q_g2_r2 q_g3_r1 q_g3_r2)) { my ($row) = read_table("$data/$f.txt"); $T{$f} = $row->[2]; }
sub Xi { my ($h, $tau, $s, @j) = @_; my $x = -$tau/(1+$tau); my $c = 1 + $x; my $r = scalar @j; my %v = (x => "$x", s => "$s"); $v{"j$_"} = "$j[$_-1]" for 1..$r;
  return $c**(3*$h - 3 + 2*$r) * eval_expr($T{"q_g${h}_r$r"}, \%v); }
print "(a) Theorem B at genus 2 from the tables:\n";
for my $ts ('2/7', '3/2', '-7/3', '5') { my $tau = R($ts); my $sig = $tau*(1+$tau); my $sp = 2*$tau + 1;
  my $Pi = sub { my $j = shift; my $v = R(1); for my $b (1..$j-1) { $v *= ($j*$tau + $b) / $b; } return $v; };
  my $x11 = -(1 + 1/$sig)/24; my $x12 = -1/(24*$sig);
  my $B3 = Xi(2, $tau, -1, 3, 2) - 4*$Pi->(4)/($Pi->(3)*$sp) * Xi(2, $tau, -1, 4) + (2*$sig/$Pi->(3)) * $x12 * (2 + $x11);
  my $B2 = Xi(2, $tau, -1, 2, 2) - 3*$Pi->(3)/$sp**2 * Xi(2, $tau, -1, 3) + $sig/(2*$sp**2) * $x11 * (1 + $x11);
  print "   tau=$ts: R=(3): $B3   R=(2): $B2   [Xi_2(3,2) = ", Xi(2, $tau, -1, 3, 2), "]\n"; }
print "(b) order of vanishing of Q_h(s;R) at s = -1 (tau = 2/7 and 5/3; derivative by exact difference quotients of the polynomial):\n";
for my $case ([2, [2]], [3, [2]], [3, [3]], [3, [4]], [3, [2,2]], [3, [3,2]]) { my ($h, $j) = @$case;
  for my $ts ('2/7', '5/3') { my $tau = R($ts); my $N = 3*$h - 3 + scalar(@$j);
    # Q is a polynomial in s of degree <= N: interpolate at s = -1, 0, 1, ..., N-1 and expand around s = -1
    my @pts = map { [ R($_), Xi($h, $tau, $_, @$j) ] } (-1..$N-1);
    # Newton divided differences at nodes -1,0,..: coefficient of (s+1) is the first divided difference structure; compute Taylor coefficients at -1
    my @nodes = map { $_->[0] } @pts; my @dd = map { $_->[1] } @pts; my $n = scalar @dd;
    for my $k (1..$n-1) { for (my $i = $n-1; $i >= $k; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($nodes[$i] - $nodes[$i-$k]); } }
    # Q(s) = dd0 + dd1 (s+1) + dd2 (s+1) s + ... ; value and first derivative at s = -1:
    my $q0 = $dd[0]; my $q1 = R(0); { my @poly = (R(1)); for my $k (1..$n-1) { # product (s - nodes[0..k-1]) evaluated derivative at -1
        my @np = map { R(0) } 0..scalar(@poly); for my $i (0..$#poly) { $np[$i+1] += $poly[$i]; $np[$i] -= $poly[$i] * $nodes[$k-1]; } @poly = @np;
        # derivative at s = -1 of poly
        my $d = R(0); for my $i (1..$#poly) { $d += $i * $poly[$i] * (R(-1) ** ($i-1)); } $q1 += $dd[$k] * $d; } }
    printf "   h=%d R=(%s) tau=%s: Q(-1) = %s, Q'(-1) = %s\n", $h, join(',', @$j), $ts, $q0, $q1; } }
print "(c) lambda_1 component of Xi_2(1/2, 5/2): sigma^2 Xi_2 = X0 sigma^2 - X1 sigma + X2; solve from three tau values\n";
{ my @rows; for my $ts ('2/7', '3/2', '5') { my $tau = R($ts); my $sig = $tau*(1+$tau); push @rows, [ $sig, $sig**2 * Xi(2, $tau, -1, R('1/2'), R('5/2')) ]; }
  # solve a sigma^2 + b sigma + c = v
  my ($s1, $v1) = @{$rows[0]}; my ($s2, $v2) = @{$rows[1]}; my ($s3, $v3) = @{$rows[2]};
  my $a = (($v1-$v2)/($s1-$s2) - ($v2-$v3)/($s2-$s3)) / ($s1 - $s3); my $b = ($v1-$v2)/($s1-$s2) - $a*($s1+$s2); my $cc = $v1 - $a*$s1**2 - $b*$s1;
  print "   Xi^(0) = $a, Xi^(1) = ", -$b, ", Xi^(2) = $cc   (lambda_g formula b_2 (J-1)(J-2)(J-3) with J = 3: 0)\n"; }
