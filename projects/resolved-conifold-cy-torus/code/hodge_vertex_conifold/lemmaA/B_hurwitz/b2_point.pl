#!/usr/bin/perl
# One point: Xi^(0)_g(R) from connected Hurwitz numbers (finite formula B2) and from the Toda recursion (if delta != 0).
# Usage: perl b2_point.pl g j1 j2 ...
use strict; use warnings; use lib '.'; use HB; use HBconn; use Rec; use Math::BigRat;
my ($g, @R) = @ARGV; my $t0 = time;
my $r = @R; my $J = 0; $J += $_ for @R; my $D = $J - $r; my $delta = $D - 2*$g + 2;
my $xh = HBconn::xi0(\@R, $g);
my $xr = $delta == 0 ? "undetermined (delta=0)" : Rec::xi($g, \@R);
my $ok = $delta == 0 ? "" : ($xh == $xr ? "ok" : "MISMATCH");
printf "g=%d R=(%s) D=%d delta=%d: Hurwitz/B2 %s | recursion %s %s [%d s]\n", $g, join(",", @R), $D, $delta, $xh, $xr, $ok, time - $t0;
