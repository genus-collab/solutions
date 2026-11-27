#!/usr/bin/perl
# The Toda leading-coefficient recursion (Rec.pm) against the exact Hodge data q_g<g>_r<r>.txt (k = 0 component at s = -1).
# Usage: perl rec_vs_data.pl
use strict; use warnings; use lib '.'; use HB; use Data; use Rec; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my ($n, $bad) = (0, 0);
for my $spec (['q_g2_r1.txt', 2, 1, 12], ['q_g2_r2.txt', 2, 2, 12], ['q_g3_r1.txt', 3, 1, 14], ['q_g3_r2.txt', 3, 2, 12],
              ['q_g4_r1.txt', 4, 1, 16], ['q_g4_r2.txt', 4, 2, 13]) {
  my ($file, $g, $r, $JMAX) = @$spec; next unless -e "../data/$file";
  my $ast = Data::load($file)->[0]{ast};
  my @Rs = $r == 1 ? (map { [$_] } 2..$JMAX) : (map { my $a = $_; map { [$a, $_] } $a..($JMAX - $a) } 2..int($JMAX / 2));
  for my $R (@Rs) { my $J = 0; $J += $_ for @$R; my $D = $J - $r; my $delta = $D - 2*$g + 2;
    my %env = (s => BR(-1)); $env{"j$_"} = BR($R->[$_ - 1]) for 1..$r;
    my $xd = Data::components($ast, $g, $r, \%env)->[0];
    if ($delta == 0) { printf "g=%d R=(%s) D=%d delta=0: data %s (not determined by the recursion)\n", $g, join(",", @$R), $D, $xd; next }
    my $xr = Rec::xi($g, $R); $n++; my $ok = $xr == $xd ? "ok" : "MISMATCH"; $bad++ if $ok ne "ok";
    printf "g=%d R=(%s) D=%d delta=%d: recursion %s, Hodge data %s %s\n", $g, join(",", @$R), $D, $delta, $xr, $xd, $ok } }
print "rec_vs_data: $n comparisons, $bad mismatches\n";
