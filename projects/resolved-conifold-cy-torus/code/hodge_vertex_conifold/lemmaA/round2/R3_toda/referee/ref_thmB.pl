#!/usr/bin/perl
# Referee check of Route A's Theorem B (cut-and-join 2-insertion relation), own implementation, own data access (ref_data.pm):
#  Xi_g(R+2) - Sum_i gam_{j_i} Xi_g(R - j_i + (j_i+1)) + Sum_{i<k} gam_{j_i j_k} Xi_g(R - j_i - j_k + (j_i+j_k))
#   + (sigma/(2 sigma')) Sum_i Sum_{al+be=j_i, al,be>=1} al be Pi_al Pi_be/Pi_{j_i} [ Xi_{g-1}(al,be,R-j_i)
#        + Sum_{g1+g2=g} Sum_{labelled I+J = R-j_i} Xi_{g1}(al,R_I) Xi_{g2}(be,R_J) ] = 0,
#  gam_j = (j+1)Pi_{j+1}/(Pi_j sigma'), gam_jk = (j+k)Pi_{j+k}/(Pi_j Pi_k sigma'), sigma' = 2tau+1.
# Usage: perl ref_thmB.pl gmin gmax rmax Jmax [taus]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use ref_data;
$| = 1;
sub R { Math::BigRat->new(@_) }
my ($gmin, $gmax, $rmax, $Jmax, $ts) = @ARGV; $ts //= "5/11,-7/3,4/9";
my @taus = map { R($_) } split /,/, $ts;
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub Pi { my ($k, $t) = @_; my $p = R(1); $p *= ($k*$t + $_) for 1..$k-1; $p / fact($k-1) }
sub X { my ($g, $T, $t) = @_; return R(0) if $g < 0; my @T = sort { $a <=> $b } @$T;
  if ($g >= 1 && !@T) { die "empty profile at genus $g" } return ref_data::xi_num($g, \@T, $t); }
sub relB { my ($g, $R, $t) = @_; my @R = @$R; my $r = @R; my $sig = $t*($t+1); my $sp = 2*$t + 1;
  my $tot = X($g, [@R, 2], $t);
  for my $i (0..$r-1) { my @c = @R; my $j = $c[$i]; $c[$i] = $j + 1; $tot -= ($j+1)*Pi($j+1,$t)/(Pi($j,$t)*$sp) * X($g, \@c, $t); }
  for my $i (0..$r-1) { for my $k ($i+1..$r-1) { my @c = map { $R[$_] } grep { $_ != $i && $_ != $k } 0..$r-1; my ($a, $b) = ($R[$i], $R[$k]);
      $tot += ($a+$b)*Pi($a+$b,$t)/(Pi($a,$t)*Pi($b,$t)*$sp) * X($g, [@c, $a+$b], $t); } }
  my $q = R(0);
  for my $i (0..$r-1) { my $j = $R[$i]; my @rem = map { $R[$_] } grep { $_ != $i } 0..$r-1; my $nr = @rem;
    for my $al (1..$j-1) { my $be = $j - $al; my $w = $al*$be*Pi($al,$t)*Pi($be,$t)/Pi($j,$t);
      my $s = ($g - 1 >= 0) ? X($g-1, [$al, $be, @rem], $t) : R(0);
      for my $g1 (0..$g) { my $g2 = $g - $g1;
        for my $mask (0..(1<<$nr)-1) { my (@I, @J); for my $m (0..$nr-1) { if ($mask & (1<<$m)) { push @I, $rem[$m] } else { push @J, $rem[$m] } }
          $s += X($g1, [$al, @I], $t) * X($g2, [$be, @J], $t); } }
      $q += $w * $s; } }
  return $tot + $sig/(2*$sp) * $q; }
my @Rs; my $gen; $gen = sub { my ($pre, $min, $J) = @_; for (my $j = $min; $J + $j <= $Jmax; $j++) { my @p = (@$pre, $j); push @Rs, [@p]; $gen->(\@p, $j, $J + $j) if @p < $rmax; } }; $gen->([], 2, 0);
my ($n, $bad) = (0, 0);
for my $g ($gmin..$gmax) { for my $R (@Rs) { my @res = map { relB($g, $R, $_) } @taus; my $ok = !grep { $_ != 0 } @res; $n++; $bad++ unless $ok;
  my $D = 1; $D += $_ - 1 for @$R; printf "g=%d R=(%s) [level of R+2: %d, excess %d] %s\n", $g, join(',', @$R), $D, 2*$g-2-$D, ($ok ? "OK" : "FAIL @res"); } }
printf "SUMMARY Theorem B (referee implementation): %d relations at %d values of tau, %d failures\n", $n, scalar(@taus), $bad;
