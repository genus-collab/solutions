#!/usr/bin/perl
# Genus-7 point values (all components) needed for the genus-7 test of the leading-order mKP relation.
# Usage: HODGE_NCH=4 perl g7_points.pl 7 "13;2,13;14;2,14"   -> prints "g | parts | v-poly coefficients"
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use Time::HiRes qw(time);
$| = 1;
my ($g, $list) = @ARGV; my $t0 = time; Hodge::setup($g); printf "# setup genus %d done %.1f s (NCH=%s)\n", $g, time - $t0, $ENV{HODGE_NCH} // 3;
for my $pt (split /;/, $list) { my @T = split /,/, $pt; my $t1 = time;
  my $v = HodgePt::WP([], \@T);
  printf "%d | %s | %s\n", $g, join(',', @T), join(' ; ', map { "$_" } @$v);
  printf "# time %.1f s\n", time - $t1; }
