#!/usr/bin/perl
# Check the Toda leading-coefficient recursion (Rec.pm) against the independent finite formula (B2) built from
# connected Hurwitz numbers (HBconn.pm), for all multisets R of parts >= 2 with J <= JMAX at genus g, delta != 0.
# Usage: perl rec_check.pl g JMAX [RMAX]
use strict; use warnings; use lib '.'; use HB; use HBconn; use Rec; use Math::BigRat;
my ($g, $JMAX, $RMAX) = @ARGV; $RMAX //= 9;
my @Rs; my $gen; $gen = sub { my ($pre, $min, $left) = @_; push @Rs, [@$pre] if @$pre;
  return if @$pre >= $RMAX; for my $j ($min..$left) { $gen->([@$pre, $j], $j, $left - $j) } };
$gen->([], 2, $JMAX);
my ($n, $bad, $zeros) = (0, 0, 0);
for my $R (sort { @$a <=> @$b || join(",", @$a) cmp join(",", @$b) } @Rs) {
  my $r = @$R; my $J = 0; $J += $_ for @$R; my $D = $J - $r; my $delta = $D - 2*$g + 2;
  next if $delta == 0;
  my $x_rec = Rec::xi($g, $R); my $x_hur = HBconn::xi0($R, $g);
  $n++; my $ok = $x_rec == $x_hur ? "ok" : "MISMATCH"; $bad++ if $ok ne "ok";
  my $range = $D <= 2*$g - 3 ? "range" : ($D == 2*$g - 1 ? "D=2g-1" : "");
  $zeros++ if $range && $x_hur == 0;
  printf "g=%d R=(%s) D=%d delta=%d %-7s recursion %s | Hurwitz %s %s\n", $g, join(",", @$R), $D, $delta, $range, $x_rec, $x_hur, $ok }
print "rec_check g=$g JMAX=$JMAX: $n comparisons, $bad mismatches\n";
