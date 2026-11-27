#!/usr/bin/perl
# ref_cj2.pl MODE G TAU JMAX : Route A's Theorem B (MV cut-and-join at the dilaton point) for two-leaf R = (j1, j2), 1 <= j2 <= j1 <= JMAX:
#   E_G(R) = Xi_G(j1,j2,2) - gam_{j1} Xi_G(j1+1,j2) - gam_{j2} Xi_G(j1,j2+1) + gam_{j1j2} Xi_G(j1+j2)
#     + (sigma/(2 sigma')) Sum_{i} Sum_{al+be=j_i} al be (Pi_al Pi_be/Pi_{j_i}) [ Xi_{G-1}(al,be,j_o) + Sum_{g1+g2=G} ( Xi_{g1}(al,j_o) Xi_{g2}(be) + Xi_{g1}(al) Xi_{g2}(be,j_o) ) ]
#   (j_o = the other part), gam_j = (j+1)Pi_{j+1}/(Pi_j sigma'), gam_{jk} = (j+k)Pi_{j+k}/(Pi_j Pi_k sigma'); genus 0: Xi_0(a) = 1/(a(a+1)), Xi_0(a,b) = 1/(a+b).
# MODE engine: all values from TR-free engine tables; MODE formula: one-leaf (6.1), two-leaf RefXY::xi2, three-leaf RefXY3::xi3 (the R2 Theorem 1).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use RefXY; use RefXY3; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($mode, $G, $tau, $jmax, $jmin) = @ARGV; $tau = Q($tau); $jmin //= 1; my $sig = $tau * (1 + $tau); my $sp = 2*$tau + 1;
my $EF = "$FindBin::Bin/../../../E_formulas"; my %T;
sub tab { my ($g, $r) = @_; my $k = "$g,$r"; return $T{$k} if exists $T{$k};
  for my $f ("$EF/engine_h${g}_r${r}.txt", "$EF/referee/ref_engine_h${g}_r${r}.txt", "$FindBin::Bin/eng/ref_engine_h${g}_r${r}.txt") {
    if (-e $f) { my $t = RefR2::load_table($f, $g, $r); if (%$t) { $T{$k} = $t; return $t; } } } die "no engine table for genus $g r=$r"; }
my %cache; my %M;
sub X { my ($g, @m) = @_; @m = sort { $b <=> $a } @m; my $r = scalar @m; my $key = "$g|@m"; return $M{$key} if exists $M{$key}; my $v;
  if ($g == 0) { my $J = 0; $J += $_ for @m; if ($r == 1) { $v = Q(1) / ($J * ($J + 1)); } elsif ($r == 2) { $v = Q(1) / $J; } else { $v = Q(1); $v *= ($J - $_) for 1..$r-3; } }
  elsif ($mode eq 'engine') { $v = RefR2::table_value(tab($g, $r), $tau, @m); }
  elsif ($r == 1) { $v = RefR2::xi61($g, $tau, $m[0], \%cache); }
  elsif ($r == 2) { $v = RefXY::xi2($tau, $g, @m); }
  else { $v = RefXY3::xi3($tau, $g, reverse @m); }   # ascending order: same value (checked), much faster
  $M{$key} = $v; return $v; }
sub gam1 { my $j = shift; return ($j + 1) * RefR2::Pi($j + 1, $tau) / (RefR2::Pi($j, $tau) * $sp); }
sub gam2 { my ($j, $k) = @_; return ($j + $k) * RefR2::Pi($j + $k, $tau) / (RefR2::Pi($j, $tau) * RefR2::Pi($k, $tau) * $sp); }
my ($n, $bad) = (0, 0);
for my $j1 ($jmin..$jmax) { for my $j2 (1..$j1) { my @R = ($j1, $j2);
  my $E = X($G, $j1, $j2, 2) - gam1($j1) * X($G, $j1 + 1, $j2) - gam1($j2) * X($G, $j1, $j2 + 1) + gam2($j1, $j2) * X($G, $j1 + $j2);
  for my $i (0, 1) { my $j = $R[$i]; my $jo = $R[1 - $i];
    for my $al (1..$j-1) { my $be = $j - $al; my $cf = $sig / (2*$sp) * $al * $be * RefR2::Pi($al, $tau) * RefR2::Pi($be, $tau) / RefR2::Pi($j, $tau);
      my $s = X($G - 1, $al, $be, $jo);
      for my $g1 (0..$G) { my $g2 = $G - $g1; $s += X($g1, $al, $jo) * X($g2, $be) + X($g1, $al) * X($g2, $be, $jo); }
      $E += $cf * $s; } }
  $n++; $bad++ unless $E->is_zero;
  printf "%s G=%d tau=%s R=(%d,%d): Xi_G(R,2)=%s E=%s %s\n", $mode, $G, $tau, $j1, $j2, X($G, $j1, $j2, 2), $E, ($E->is_zero ? "OK" : "FAIL"); } }
printf "SUMMARY %s G=%d tau=%s j<=%d: %d relations, %d failures\n", $mode, $G, $tau, $jmax, $n, $bad;
