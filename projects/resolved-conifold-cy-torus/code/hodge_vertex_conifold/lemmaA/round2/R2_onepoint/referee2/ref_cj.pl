#!/usr/bin/perl
# ref_cj.pl MODE G TAU JMIN JMAX : TR-free consistency test with Route A's Theorem B (Marino-Vafa cut-and-join at the dilaton point),
# one-leaf case R = (j):
#   E_G(j) = Xi_G(j,2) - gam_j Xi_G(j+1)
#            + (sigma/(2 sigma')) Sum_{al+be=j} al be (Pi_al Pi_be/Pi_j) [ Xi_{G-1}(al,be) + Sum_{g1=0}^{G} Xi_{g1}(al) Xi_{G-g1}(be) ]  = 0,
#   gam_j = (j+1) Pi_{j+1}/(Pi_j sigma'), sigma' = 2 tau + 1;  genus-0 values Xi_0(a) = 1/(a(a+1)).
# MODE = engine  : every Xi from the TR-free engine tables (validates the relation and its conventions);
# MODE = formula : one-leaf values from formula (6.1), two-leaf values from Theorem 1 (n = 2, RefXY::xi2) -- the R2 predictions.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use RefXY; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($mode, $G, $tau, $jmin, $jmax) = @ARGV; $tau = Q($tau); my $sig = $tau * (1 + $tau); my $sp = 2*$tau + 1;
my $EF = "$FindBin::Bin/../../../E_formulas"; my %T;
sub tab { my ($g, $r) = @_; my $k = "$g,$r"; return $T{$k} if exists $T{$k};
  my @c = ("$EF/engine_h${g}_r${r}.txt", "$EF/referee/ref_engine_h${g}_r${r}.txt", "$FindBin::Bin/eng/ref_engine_h${g}_r${r}.txt");
  for my $f (@c) { if (-e $f) { my $t = RefR2::load_table($f, $g, $r); if (%$t) { $T{$k} = $t; return $t; } } } die "no engine table for genus $g r=$r"; }
my %cache;
sub X1 { my ($g, $a) = @_; return Q(1) / ($a * ($a + 1)) if $g == 0;
  return $mode eq 'engine' ? RefR2::table_value(tab($g, 1), $tau, $a) : RefR2::xi61($g, $tau, $a, \%cache); }
sub X2 { my ($g, $a, $b) = @_; die "genus 0 two-leaf not needed" if $g == 0;
  return $mode eq 'engine' ? RefR2::table_value(tab($g, 2), $tau, $a, $b) : RefXY::xi2($tau, $g, ($a >= $b ? ($a, $b) : ($b, $a))); }
my ($n, $bad) = (0, 0);
for my $j ($jmin..$jmax) {
  my $gam = ($j + 1) * RefR2::Pi($j + 1, $tau) / (RefR2::Pi($j, $tau) * $sp);
  my $t1 = X2($G, $j, 2); my $t2 = $gam * X1($G, $j + 1);
  my $t3 = Q(0);
  for my $al (1..$j-1) { my $be = $j - $al; my $cf = $sig / (2*$sp) * $al * $be * RefR2::Pi($al, $tau) * RefR2::Pi($be, $tau) / RefR2::Pi($j, $tau);
    my $s = X2($G - 1, $al, $be); for my $g1 (0..$G) { $s += X1($g1, $al) * X1($G - $g1, $be); } $t3 += $cf * $s; }
  my $E = $t1 - $t2 + $t3; $n++;
  my $ok = $E->is_zero ? "OK" : "FAIL"; $bad++ unless $E->is_zero;
  printf "%s G=%d tau=%s j=%d: Xi_G(j,2)=%s  E=%s %s\n", $mode, $G, $tau, $j, $t1, $E, $ok; }
printf "SUMMARY %s G=%d tau=%s j=%d..%d: %d relations, %d failures\n", $mode, $G, $tau, $jmin, $jmax, $n, $bad;
