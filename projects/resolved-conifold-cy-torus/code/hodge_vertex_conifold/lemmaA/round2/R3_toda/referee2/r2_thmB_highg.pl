#!/usr/bin/perl
# Second referee: Route A's Theorem B (2-insertion relation from the MV cut-and-join at the dilaton point), used in proof (i),
# tested at genus G >= 7 for one-part R = (j) outside the range (non-trivial), mod p, with the referee's own MV engine (R2mv.pm)
# for genus >= 7 and Route E's tables for genus <= 6 (re-coded from the statement in A_mv/report.md s.3.3):
#  Xi_g(R+2) - sum_i gam_{j_i} Xi_g(R_i^+) + sum_{i<k} gam_{j_i j_k} Xi_g(R_ik)
#   + (sigma/(2 sigma')) sum_i sum_{al+be=j_i} al be (Pi_al Pi_be/Pi_{j_i}) [Xi_{g-1}(al,be,R-j_i) + sum_{g1+g2=g} sum_{I+J} Xi_g1(al,R_I) Xi_g2(be,R_J)] = 0
#  gam_j = (j+1)Pi_{j+1}/(Pi_j sigma'), sigma' = 2tau+1; unit parts by string Xi_g(1,S) = J_S Xi_g(S); Xi_g(1) = 0 (g >= 2),
#  Xi_1(1) = -(1+1/sigma)/24; genus 0: Gamma(J)/Gamma(J-r+3) for all r >= 1.
# Usage: perl r2_thmB_highg.pl prime seed G "j1;j2;..."
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use R2mv; use R2tab;
$| = 1;
my ($prime, $seed, $G, $jl) = @ARGV; R2mv::setp($prime); my $p = $prime; srand($seed); my $tau = 2 + int(rand($p - 10));
my @js = split /;/, $jl; my $jmax = 0; for (@js) { $jmax = $_ if $_ > $jmax; }
my $sig = $tau * R2mv::md($tau + 1) % $p; my $vv = R2mv::md(-R2mv::inv($sig));
# genus-G values needed: (2,j), (j+1), (a) for 2 <= a <= jmax-1 (one-leaf join factors), (a,b) with a+b = j at genus G-1 (if G-1 > 6)
my %Bs; for my $j (@js) { $Bs{"2,$j"} = [2, $j]; $Bs{$j + 1} = [$j + 1]; } for my $a (2..$jmax - 1) { $Bs{$a} = [$a]; }
if ($G - 1 > 6) { for my $j (@js) { for my $al (2..$j - 2) { my $be = $j - $al; next if $be < $al; $Bs{"$al,$be"} = [$al, $be]; } } }
my @B = values %Bs; my $t0 = time; my $res = R2mv::compute($G, \@B, $tau, report => 1);
my %XG; for my $k (keys %Bs) { my @U = sort { $a <=> $b } @{$Bs{$k}}; my $kd = join(',', sort { $b <=> $a } @U); for my $g (1..$G) { $XG{"$g|" . join(',', @U)} = $res->{$kd}{Xi}{$g} if defined $res->{$kd}{Xi}{$g}; } }
my $ncross = 0; my $badcross = 0;
sub tabval { my ($g, $U) = @_; my $c = R2tab::comps($g, $U); my $s = 0; my $q = 1; for my $k (0..$#$c) { $s = ($s + R2mv::rat($c->[$k]) * $q) % $p; $q = $q * $vv % $p; } $s }
sub X { my ($g, $T) = @_; my @T = sort { $a <=> $b } @$T; return 0 if $g < 0;
  if ($g == 0) { return @T ? R2mv::rat(R2tab::xi0(\@T)) : die "X_0(empty)"; }
  if (@T && $T[0] == 1) { shift @T; my $Js = 0; $Js += $_ for @T; if (!@T) { return $g == 1 ? R2mv::md(-(1 + R2mv::inv($sig)) * R2mv::inv(24)) : 0; } return $Js * X($g, \@T) % $p; }
  my $key = "$g|" . join(',', @T);
  if ($g <= 6) { my $t = tabval($g, \@T); if (exists $XG{$key}) { $ncross++; $badcross++ if $XG{$key} != $t; } return $t; }
  die "missing genus-$g value (@T)" unless exists $XG{$key}; $XG{$key} }
sub Pi { R2mv::Pi($_[0], $tau) }
my $sp = R2mv::md(2*$tau + 1);
for my $j (@js) { my $tot = X($G, [2, $j]);
  $tot = R2mv::md($tot - ($j + 1) * Pi($j + 1) % $p * R2mv::inv(Pi($j) * $sp % $p) % $p * X($G, [$j + 1]));
  my $q = 0;
  for my $al (1..$j - 1) { my $be = $j - $al; my $w = $al * $be % $p * Pi($al) % $p * Pi($be) % $p * R2mv::inv(Pi($j)) % $p;
    my $s = X($G - 1, [$al, $be]);
    for my $g1 (0..$G) { $s = ($s + X($g1, [$al]) * X($G - $g1, [$be])) % $p; }
    $q = ($q + $w * $s) % $p; }
  $tot = ($tot + $sig * R2mv::inv(2 * $sp % $p) % $p * $q) % $p;
  my $D = $j; printf "Theorem B, g=%d R=(%d) [R+2 = (2,%d): D=%d, excess %d]: residual %d %s\n", $G, $j, $j, $D, 2*$G - 2 - $D, $tot, ($tot ? "FAIL" : "OK"); }
printf "# genus <= 6 engine values cross-checked against tables: %d, mismatches %d; total %d s\n", $ncross, $badcross, time - $t0;
