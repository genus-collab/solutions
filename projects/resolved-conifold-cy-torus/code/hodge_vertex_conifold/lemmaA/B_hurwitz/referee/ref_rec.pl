#!/usr/bin/perl
# Referee: (1) the referee's (REC) implementation versus the author's Rec.pm on a large range (implementation cross-check);
# (2) zeros predicted by Theorem A inside the recursion; (3) Theorem B, both closed forms, versus (REC) for r = 1;
# (4) Lemma 6.3 (pure-unit cancellation u_m = 0) with L_{g,k,empty} = (2g-1) B_2g (2g+2k-3)!/(2g)!.
# Usage: perl ref_rec.pl GMAX JMAX RMAX
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use lib "$FindBin::Bin/.."; use RefRec; use Rec; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my ($GMAX, $JMAX, $RMAX) = @ARGV; $GMAX //= 5; $JMAX //= 12; $RMAX //= 6;
my ($n, $bad, $z, $zbad) = (0, 0, 0, 0);
my @Rs; my $gen; $gen = sub { my ($pre, $min, $left) = @_; push @Rs, [@$pre] if @$pre; return if @$pre >= $RMAX; for my $j ($min..$left) { $gen->([@$pre, $j], $j, $left - $j) } };
$gen->([], 2, $JMAX);
for my $g (1..$GMAX) { for my $R (sort { @$a <=> @$b || join(",", @$a) cmp join(",", @$b) } @Rs) {
  my $r = @$R; my $J = 0; $J += $_ for @$R; my $D = $J - $r; my $delta = $D - 2 * $g + 2; next if $delta == 0;
  my $mine = RefRec::xi($g, @$R); my $his = Rec::xi($g, [@$R]); $n++;
  if (!defined $his || $mine != $his) { $bad++; print "MISMATCH g=$g R=(@$R): referee $mine, author ", $his // 'undef', "\n" }
  if ($D <= 2 * $g - 3 || $D == 2 * $g - 1) { $z++; if (!$mine->is_zero) { $zbad++; print "THEOREM A VIOLATED in REC: g=$g R=(@$R) value $mine\n" } } } }
print "(1) referee REC vs author REC (g <= $GMAX, J <= $JMAX, r <= $RMAX, delta != 0): $n comparisons, $bad mismatches\n";
print "(2) Theorem-A points (D <= 2g-3 or D = 2g-1) among them: $z, non-zero: $zbad\n";
my ($m3, $b3) = (0, 0);
for my $g (1..10) { for my $j (2..3 * $g + 2) { my $f1 = RefRec::thmB_form1($g, $j); my $f2 = RefRec::thmB_form2($g, $j);
  $m3++; if ($f1 != $f2) { $b3++; print "THM B forms differ g=$g j=$j: $f1 vs $f2\n" }
  my $delta = $j - 2 * $g + 1; next if $delta == 0;
  my $rec = RefRec::xi($g, $j); $m3++; if ($rec != $f1) { $b3++; print "THM B vs REC g=$g j=$j: $f1 vs $rec\n" }
  if ($j == 2 * $g && !$f1->is_zero) { $b3++; print "extra zero j=2g fails g=$g\n" } } }
print "(3) Theorem B (two closed forms) vs REC, g <= 10: $m3 comparisons, $b3 failures\n";
my ($m4, $b4) = (0, 0);
for my $m (1..25) { my $u = BR(0); for my $g (0..$m) { my $k = $m + 1 - $g; next if $g == 0 && $k == 1; $u += BR(2) / BR(RefRec::fac(2 * $k)) * RefRec::L_empty($g, $k) }
  $m4++; if (!$u->is_zero) { $b4++; print "u_$m = $u (non-zero)\n" } }
# consistency of L_{g,k,empty} with D^{2k} applied to chi(M_g) T^{-(2g-2)}: (-(2g-2))_(2k) * B_2g/(2g(2g-2))
for my $g (2..8) { for my $k (1..5) { my $chi = RefRec::bern(2 * $g) / BR(2 * $g * (2 * $g - 2)); my $L = RefRec::ff(-(2 * $g - 2), 2 * $k) * $chi;
  $m4++; if ($L != RefRec::L_empty($g, $k)) { $b4++; print "L formula mismatch g=$g k=$k\n" } } }
print "(4) Lemma 6.3 u_m = 0 (m <= 25) and L-formula consistency: $m4 checks, $b4 failures\n";
