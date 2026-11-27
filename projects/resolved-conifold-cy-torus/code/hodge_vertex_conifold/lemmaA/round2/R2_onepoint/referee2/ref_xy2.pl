#!/usr/bin/perl
# ref_xy2.pl TAU G JMAX [TABLE] : Theorem 1 (n = 2) by RefXY::xi2 (referee implementation) for 1 <= m2 <= m1 <= JMAX at integer TAU,
# compared with a TR-free engine table (r = 2) when given; also checks the symmetry m1 <-> m2 of the iterated residue and the
# vanishing range |mu| <= 2G-1.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use RefXY; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($tau, $g, $jmax, $file) = @ARGV; $tau = Q($tau);
my $E = defined $file ? RefR2::load_table($file, $g, 2) : undef;
my ($n, $bad, $asym, $rng, $nz) = (0, 0, 0, 0, 0);
for my $m1 (1..$jmax) { for my $m2 (1..$m1) {
  my $x = RefXY::xi2($tau, $g, $m1, $m2); my $y = RefXY::xi2($tau, $g, $m2, $m1); $asym++ if $x != $y;
  $n++; $nz++ unless $x->is_zero;
  if ($m1 + $m2 <= 2*$g - 1 && !$x->is_zero) { $rng++; print "RANGE VIOLATION ($m1,$m2): $x\n"; }
  if (defined $E) { my $e = RefR2::table_value($E, $tau, $m1, $m2); if ($e != $x) { $bad++; print "MISMATCH tau=$tau g=$g ($m1,$m2): xy $x engine $e\n" if $bad <= 10; } } } }
printf "ref_xy2 tau=%s g=%d m<=%d: %d pairs (%d non-zero), asymmetric %d, range violations %d, %s\n", $tau, $g, $jmax, $n, $nz, $asym, $rng,
  (defined $E ? "engine mismatches $bad" : "no engine table");
