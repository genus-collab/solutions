#!/usr/bin/perl
# ref_xy3.pl TAU G JMAX [TABLE] [PERMCHECK] : Theorem 1 at n = 3 (RefXY3, referee implementation) for 1 <= m3 <= m2 <= m1 <= JMAX at integer
# TAU, compared with a TR-free engine table (r = 3); range check |mu| <= 2G; optional symmetry check under a permutation of (m1,m2,m3).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use RefXY; use RefXY3; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($tau, $g, $jmax, $file, $perm) = @ARGV; $tau = Q($tau);
my $E = (defined $file && $file ne '-') ? RefR2::load_table($file, $g, 3) : undef;
my ($n, $bad, $asym, $rng, $nz) = (0, 0, 0, 0, 0);
for my $m1 (1..$jmax) { for my $m2 (1..$m1) { for my $m3 (1..$m2) {
  my $x = RefXY3::xi3($tau, $g, $m1, $m2, $m3); $n++; $nz++ unless $x->is_zero;
  if ($perm) { my $y = RefXY3::xi3($tau, $g, $m3, $m1, $m2); my $z = RefXY3::xi3($tau, $g, $m2, $m3, $m1); $asym++ if $x != $y || $x != $z; }
  if ($m1 + $m2 + $m3 <= 2*$g && !$x->is_zero) { $rng++; print "RANGE VIOLATION ($m1,$m2,$m3): $x\n"; }
  if (defined $E) { my $e = RefR2::table_value($E, $tau, $m1, $m2, $m3); if ($e != $x) { $bad++; print "MISMATCH tau=$tau g=$g ($m1,$m2,$m3): xy $x engine $e\n" if $bad <= 10; } } } } }
printf "ref_xy3 tau=%s g=%d m<=%d: %d triples (%d non-zero), %s, range violations %d, %s\n", $tau, $g, $jmax, $n, $nz,
  ($perm ? "asymmetric $asym" : "no symmetry check"), $rng, (defined $E ? "engine mismatches $bad" : "no engine table");
