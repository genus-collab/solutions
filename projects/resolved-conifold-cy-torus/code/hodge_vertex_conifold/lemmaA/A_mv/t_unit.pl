#!/usr/bin/perl
# t_unit.pl -- consistency: unit-leg generating function (mv_unit, p_1 = zeta x) versus the plain cumulant formula (mv_conn)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my $tau = R('3/2'); my $sig = $tau*(1+$tau); my $bad = 0; my $n = 0;
for my $Rp ([2], [3], [2,2], [3,2]) { my $gmax = 2; my $A = 3;
  my $M = mv_unit($Rp, $tau, $gmax, $A);
  for my $a (0..$A) { my @mu = (@$Rp, (1) x $a); my $sz = 0; $sz += $_ for @mu; my $l = scalar @mu;
    my $P = 2*$gmax - 2 + $l + $sz; my $s = mv_conn(\@mu, $tau, $P);
    my $Pi = R(1); for my $j (@mu) { for my $b (1..$j-1) { $Pi *= ($j*$tau + $b); } $Pi /= fact($j-1); }
    for my $g (0..$gmax) { my $idx = 2*$g - 2 + $l + $sz; next if $idx < 0 || !defined $M->[$g][$a];
      my $Ig = R((-1)**$g) * $s->[$idx] * mvlib::autf(@mu) / ($sig ** ($l-1) * $Pi); my $Mg = $Ig / ((-$sig)**$g);
      $n++; if ($Mg != $M->[$g][$a]) { $bad++; print "MISMATCH R=(@$Rp) a=$a g=$g: unit=$M->[$g][$a] conn=$Mg\n"; } } } }
print "mv_unit vs mv_conn: $n comparisons, $bad mismatches\n";
