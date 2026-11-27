#!/usr/bin/perl
# ref_cmp61_engine.pl TABLE h : componentwise EXACT comparison of formula (6.1) with a TR-free engine table (r = 1).
#   (6.1) gives sigma^h Xi_h(j) = Sum_k (-1)^k sigma^{h-k} Xi^(k)(j) as a polynomial in (sigma, j); the engine gives each Xi^(k)(j) as a
#   polynomial in j.  Every coefficient (k, j^b) is compared, so a match is an identity in (tau, j).
# Also prints the gap value sigma^h F_h, F_h = Xi_h(0)/(2-2h), and the zero range of the one-leaf polynomial.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($file, $h) = @ARGV;
my $E = RefR2::load_table($file, $h, 1);
my $F = RefR2::comps61($h);
my ($n, $bad) = (0, 0);
my %keys; for my $k (keys %$E) { $keys{"$k|$_"} = 1 for keys %{$E->{$k}}; } for my $k (keys %$F) { $keys{"$k|$_"} = 1 for keys %{$F->{$k}}; }
for my $kk (sort keys %keys) { my ($k, $b) = split /\|/, $kk; my $e = $E->{$k}{$b} // Q(0); my $f = $F->{$k}{$b} // Q(0); $n++;
  if ($e != $f) { $bad++; print "  MISMATCH k=$k j^$b: engine $e  (6.1) $f\n" if $bad <= 20; } }
printf "genus %d: %d (component, j-power) coefficients compared, %d mismatches; engine components present: %s\n", $h, $n, $bad, join(',', sort { $a <=> $b } keys %$E);
# zero range and gap at a few framings, straight from the engine table
my @taus = map { Q($_) } qw(1 2 3 1/3 -3/7 5/2);
for my $tau (@taus) { my @z; for my $j (1..2*$h) { my $v = RefR2::table_value($E, $tau, $j); push @z, $j if $v->is_zero; }
  my $x0 = RefR2::table_value($E, $tau, 0); my $sig = $tau*(1+$tau); my $Fh = $x0 / (2 - 2*$h);
  # gap formula: sigma^h F_h = -(2h-3)! [lam^{2h}] S(e1 lam) S(e2 lam)  (j = 0 of the f-series)
  my @f = RefR2::fseries($h); my $g = -RefR2::fact(2*$h-3) * RefR2::eval_poly($f[$h], $sig, Q(0));
  printf "  tau=%s: engine zeros at j = %s ; sigma^h F_h engine %s, gap formula %s %s\n", $tau, join(',', @z), $sig**$h * $Fh, $g, ($sig**$h * $Fh == $g ? "OK" : "MISMATCH"); }
