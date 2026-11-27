#!/usr/bin/perl
# ref_string.pl TAU G JMAX : string-equation test of Theorem 1 across n:  Xi_G(m, 1) = m Xi_G(m)  (two-point formula vs (6.1)),
# and Xi_G(m1, m2, 1) = (m1 + m2) Xi_G(m1, m2) (three-point vs two-point) for small parts.  Both sides from the R2 formulas.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use RefXY; use RefXY3; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($tau, $G, $jmax, $j3) = @ARGV; $tau = Q($tau); $j3 //= 0; my %c; my ($n, $bad) = (0, 0);
for my $m (1..$jmax) { my $a = RefXY::xi2($tau, $G, $m, 1); my $b = $m * RefR2::xi61($G, $tau, $m, \%c); $n++; if ($a != $b) { $bad++; print "FAIL (m,1) m=$m: $a vs $b\n"; } }
for my $m1 (1..$j3) { for my $m2 (1..$m1) { my $a = RefXY3::xi3($tau, $G, $m1, $m2, 1); my $b = ($m1 + $m2) * RefXY::xi2($tau, $G, $m1, $m2); $n++; if ($a != $b) { $bad++; print "FAIL (m1,m2,1)=($m1,$m2): $a vs $b\n"; } } }
print "string test tau=$tau G=$G: $n identities, $bad failures\n";
