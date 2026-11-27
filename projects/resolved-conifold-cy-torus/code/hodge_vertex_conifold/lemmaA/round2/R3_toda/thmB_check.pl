#!/usr/bin/perl
# Exact check of Route A's Theorem B (2-insertion relation from the MV cut-and-join at the dilaton point):
#  Xi_g(R+2) - Sum_i gam_{j_i} Xi_g(R_i^+) + Sum_{i<k} gam_{j_i j_k} Xi_g(R_ik)
#   + (sigma/(2 sigma')) Sum_i Sum_{al+be=j_i} al be (Pi_al Pi_be/Pi_{j_i}) [ Xi_{g-1}(al,be,R-j_i)
#        + Sum_{g1+g2=g} Sum_{I+J=R-j_i} Xi_{g1}(al,R_I) Xi_{g2}(be,R_J) ] = 0,
#  gam_j = (j+1)Pi_{j+1}/(Pi_j sigma'), gam_jk = (j+k)Pi_{j+k}/(Pi_j Pi_k sigma'), sigma' = 2 tau + 1.
#  Unit parts by string: Xi_g(1,S) = J_S Xi_g(S); Xi_g(1) = 0 (g >= 2), Xi_1(1) = -(1+1/sigma)/24, Xi_0 for all r >= 1.
# Usage: perl thmB_check.pl gmin gmax rmax Jmax [taus]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use XiData;
$| = 1;
sub R { Math::BigRat->new(@_) }
my ($gmin, $gmax, $rmax, $Jmax, $taustr) = @ARGV; $taustr //= "2/7,3/2,-5/3";
my @taus = map { R($_) } split /,/, $taustr;
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub Pi { my ($k, $tau) = @_; my $p = R(1); for my $a (1..$k-1) { $p *= ($k*$tau + $a); } $p / fact($k-1) }
sub X { my ($g, $T, $tau) = @_;   # Xi_g(T), T may contain unit parts; returns undef for unsupported
  my @T = sort { $a <=> $b } @$T; my $sig = $tau*($tau+1);
  return R(0) if $g < 0;
  if ($g == 0) { return @T ? XiData::xi0(\@T) : undef; }
  my @ones = grep { $_ == 1 } @T; my @rest = grep { $_ >= 2 } @T;
  if (@ones) { # string equation, peel one unit part
    shift @T; my $Js = 0; $Js += $_ for @T;
    if (!@T) { return $g == 1 ? -(1 + 1/$sig)/24 : R(0); }
    return $Js * X($g, \@T, $tau); }
  return undef unless @rest;
  return XiData::xi_at($g, \@T, $tau); }
sub subsets { my @T = @_; my @out; my $r = scalar @T;
  for my $mask (0..(1<<$r)-1) { my (@a, @b); for my $i (0..$r-1) { if ($mask & (1<<$i)) { push @a, $T[$i] } else { push @b, $T[$i] } } push @out, [\@a, \@b]; } @out }
sub relB { my ($g, $Rr, $tau) = @_; my @Rr = @$Rr; my $r = scalar @Rr; my $sig = $tau*($tau+1); my $sp = 2*$tau + 1;
  my $tot = X($g, [@Rr, 2], $tau);
  for my $i (0..$r-1) { my @c = @Rr; my $j = $c[$i]; $c[$i]++; $tot -= ($j+1)*Pi($j+1,$tau)/(Pi($j,$tau)*$sp) * X($g, \@c, $tau); }
  for my $i (0..$r-1) { for my $k ($i+1..$r-1) { my @c = map { $Rr[$_] } grep { $_ != $i && $_ != $k } 0..$r-1; my ($j1, $j2) = ($Rr[$i], $Rr[$k]);
      $tot += ($j1+$j2)*Pi($j1+$j2,$tau)/(Pi($j1,$tau)*Pi($j2,$tau)*$sp) * X($g, [@c, $j1+$j2], $tau); } }
  my $q = R(0);
  for my $i (0..$r-1) { my $j = $Rr[$i]; my @rem = map { $Rr[$_] } grep { $_ != $i } 0..$r-1;
    for my $al (1..$j-1) { my $be = $j - $al; my $w = $al*$be*Pi($al,$tau)*Pi($be,$tau)/Pi($j,$tau);
      my $s = X($g-1, [$al, $be, @rem], $tau);
      for my $g1 (0..$g) { my $g2 = $g - $g1; for my $ss (subsets(@rem)) { my ($I, $J) = @$ss;
          my $a1 = X($g1, [$al, @$I], $tau); my $a2 = X($g2, [$be, @$J], $tau); $s += $a1 * $a2; } }
      $q += $w * $s; } }
  $tot += $sig/(2*$sp) * $q;
  return $tot; }
my @Rs; my $rec; $rec = sub { my ($pre, $min, $J) = @_; for (my $j = $min; $J + $j <= $Jmax; $j++) { my @p = (@$pre, $j); push @Rs, [@p]; $rec->(\@p, $j, $J + $j) if scalar(@p) < $rmax; } };
$rec->([], 2, 0);
my ($n, $bad) = (0, 0);
for my $g ($gmin..$gmax) { for my $Rr (@Rs) { my $D = 0; $D += $_ - 1 for @$Rr; $D += 1; my $exc = 2*$g - 2 - $D;
  my @res = map { relB($g, $Rr, $_) } @taus; my $ok = !(grep { $_ != 0 } @res); $n++; $bad++ unless $ok;
  printf "g=%d R={%s} [R+2: D=%d excess=%d] %s\n", $g, join(',', @$Rr), $D, $exc, ($ok ? "OK" : "FAIL @res"); } }
printf "SUMMARY Theorem B: %d relations, %d failures\n", $n, $bad;
