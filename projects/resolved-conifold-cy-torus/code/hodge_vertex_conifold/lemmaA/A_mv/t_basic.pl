#!/usr/bin/perl
# t_basic.pl -- sanity tests of mvlib.pm: characters, WK numbers, one-part MV formula.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
# characters: column orthogonality sum_nu chi_nu(mu)^2 = z_mu
for my $n (1..7) { for my $mu (partitions($n)) { my $s = 0; for my $nu (partitions($n)) { $s += chi($nu, $mu)**2; }
    die "orthogonality fails n=$n mu=@$mu" unless R($s) == zee(@$mu); } }
print "characters: column orthogonality OK for n <= 7\n";
# dimension via hook formula
for my $n (1..8) { for my $nu (partitions($n)) { my @h = hooks($nu); my $H = R(1); $H *= $_ for @h;
    die "hook formula fails @$nu" unless fact($n) / $H == chi($nu, [(1) x $n]); } }
print "hook length formula = chi(1^n): OK n <= 8\n";
print "wk checks: <tau_4>_2 = ", wk(2,4), ", <tau_1^3>_1 = ", wk(1,1,1,1), ", <tau_0 tau_2 tau_3>_2 = ", wk(2,0,2,3), " (= <tau_1 tau_3>+<tau_2 tau_2> = 3/1152*... )\n";
print "  <tau_2 tau_6>_3 = ", wk(3,2,6), " vs (5+12g(g-1))/(5 24^g g!) = ", R(5+12*3*2)/(5*R(24)**3*6), "\n";
# one-part MV: sum_g (-1)^g I_g((d)) x^{2g} = (1/d^2) S(dx) prod_{a=1}^{d-1} S(ax)/S((d tau+a)x)
my $tau = R('2/7'); my $P = 12;
for my $d (1..4) { my $s = mvlib::mv_conn([$d], $tau, $P);
  # [p_d] log R = x^{-d} s ; [x^{2g-1}] -> index 2g-1+d ;  = (-1)^g Pi_d I_g
  my $Pi = R(1); for my $b (1..$d-1) { $Pi *= ($d*$tau + $b); } $Pi /= fact($d-1);
  # closed form series of (1/d^2) S(dx) prod S(ax)/S((d tau + a) x), S(y) = exp(-sum B2k y^2k/(2k(2k)!))
  my @e = map { R(0) } 0..$P;
  for (my $k=1; 2*$k <= $P; $k++) { my $c = -bern(2*$k)/(2*$k*fact(2*$k)); my $sum = R($d)**(2*$k);
     for my $a (1..$d-1) { $sum += R($a)**(2*$k) - ($d*$tau+$a)**(2*$k); } $e[2*$k] = $c*$sum; }
  my $f = ser_exp(\@e, $P);
  my @cmp; for (my $g=0; 2*$g-1+$d <= $P && 2*$g <= $P; $g++) { my $Ig = R((-1)**$g) * $s->[2*$g-1+$d] / $Pi; my $cf = $f->[2*$g]/($d*$d) * R((-1)**$g);
     push @cmp, ($Ig == $cf ? "g$g ok" : "g$g MISMATCH $Ig vs $cf"); }
  print "d=$d: @cmp\n"; }
