#!/usr/bin/perl
# Mechanical check of the vanishing argument (Corollaries I_2 and I_P) for the leading-order relation of a Hirota member P:
# for every G in [Gmin, Gmax], every level D in [1, 2G-3] and every multiset R (parts >= 2) with D(R) = D - k + 1
# (k = weight of P; target profile R + k at level D), enumerate ALL terms of the leading-order expansion
#   (monomials of P, set partitions into blocks, sub-multiset assignments of R to blocks, (g_B, k_B) with the budget
#    Sum_B (2 g_B - 2 + |B| + k_B) = n = 2G)
# and classify each term: 'target' (linear D_{t_k}, one block, g = G, k_B = 1, T_B = R), or killed by some block:
#   ff     : U nonempty, 0 <= delta_g(U) < k_B (falling factorial vanishes)
#   H1     : U nonempty, 2 <= g < G, delta_g(U) <= -1 (Lemma A in lower genus)
#   H2     : U nonempty, g = G, D(U) < D (lower level at genus G)
#   empty0 : U empty with vanishing leading coefficient (g >= 2; g = 0, k >= 3; g = 1, k >= 1)
#   coef0  : the monomial coefficient c_alpha times combinatorial factors vanishes (never happens here)
# Any other term is a SURVIVOR (a gap). Also checks that the target occurs with a non-zero numerical coefficient.
# Usage: perl ref_count.pl EQ Gmin Gmax      (EQ = M (= D1^2 + D2), P3, P4, P5, JM2)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use ref_lib;
$| = 1;
sub R { Math::BigRat->new(@_) }
my ($EQ, $Gmin, $Gmax) = @ARGV; $EQ //= 'M';
my %EQS = (M => [[R(1), [1,1]], [R(1), [2]]], JM2 => [[R(1), [1,1,1]], [R(3), [1,2]], [R(-4), [3]]]);
for my $k (2..6) { my @P = ([R(1), [$k]], ref_lib::pk_minus_Dtilde($k)); my %h; for my $t (@P) { my $key = join(',', sort { $a <=> $b } @{$t->[1]}); $h{$key} = ($h{$key} // R(0)) + $t->[0]; }
  $EQS{"P$k"} = [ map { [$h{$_}, [split /,/, $_]] } grep { $h{$_} != 0 } keys %h ]; }
my $P = $EQS{$EQ} or die; my $wt = 0; $wt += $_ for @{$P->[0][1]};
my ($lin) = grep { @{$_->[1]} == 1 } @$P; my $clin = $lin ? $lin->[0] : R(0);
print "equation $EQ, weight $wt, coefficient of D_{t_$wt}: $clin\n";
sub setparts { my ($n) = @_; return ([]) if $n == 0; my @out;
  for my $sp (setparts($n-1)) { for my $i (0..$#$sp) { my @c = map { [@$_] } @$sp; push @{$c[$i]}, $n-1; push @out, \@c; } my @c = map { [@$_] } @$sp; push @c, [$n-1]; push @out, \@c; }
  return @out; }
# ordered decompositions of a multiset (sorted list) into q sub-multisets
sub decomp { my ($R, $q) = @_; return ([[@$R]]) if $q == 1; my %c; $c{$_}++ for @$R; my @vals = sort { $a <=> $b } keys %c; my @subs = ([]);
  for my $v (@vals) { my @n; for my $s (@subs) { push @n, [@$s, ($v) x $_] for 0..$c{$v}; } @subs = @n; }
  my @out; for my $s (@subs) { my %cs; $cs{$_}++ for @$s; my @rest; for my $v (@vals) { push @rest, ($v) x ($c{$v} - ($cs{$v} // 0)); }
    for my $d (decomp(\@rest, $q - 1)) { push @out, [$s, @$d]; } } return @out; }
sub kill_reason { my ($B, $g, $k, $T, $G, $D) = @_; my @U = (@$T, grep { $_ >= 2 } @$B);
  if (!@U) { return 'empty0' if $g >= 2 || ($g == 0 && $k >= 3) || ($g == 1 && $k >= 1); return ''; }
  my $J = 0; $J += $_ for @U; my $DU = $J - @U; my $del = $DU - 2*$g + 2;
  return 'ff' if $del >= 0 && $del < $k;
  return 'H1' if $g >= 2 && $g < $G && $del <= -1;
  return 'H2' if $g == $G && $DU < $D;
  return ''; }
my %tot; my $nsurv = 0; my $nR = 0; my $badtarget = 0;
for my $G ($Gmin..$Gmax) { for my $D (1..2*$G-3) { my $DR = $D - $wt + 1; next if $DR < 0;
  # all multisets R with parts >= 2 and D(R) = DR  <->  partitions of DR (parts j-1)
  for my $lam (ref_lib::partitions($DR)) { my @R = sort { $a <=> $b } map { $_ + 1 } @$lam; $nR++; my $n = 2*$G; my $ntarget = 0;
    for my $t (@$P) { my ($c, $idx) = @$t; my $m = @$idx;
      for my $sp (setparts($m)) { my @blocks = map { [ map { $idx->[$_] } @$_ ] } @$sp; my $q = @blocks;
        for my $dc (decomp(\@R, $q)) {
          my $rec; $rec = sub { my ($b, $left, $reasons, $gk) = @_;
            if ($b == $q) { return unless $left == 0;
              my @why = grep { $_ ne '' } @$reasons;
              my $istarget = ($m == 1 && $q == 1 && $gk->[0][0] == $G && $gk->[0][1] == 1 && scalar(@{$dc->[0]}) == scalar(@R));
              if ($istarget) { $ntarget++; $tot{target}++; return; }
              if (@why) { $tot{$why[0]}++; return; }
              $nsurv++; $tot{SURVIVOR}++;
              print "SURVIVOR: G=$G D=$D R=(@R) monomial [@$idx] blocks " . join(' ', map { "[" . join(',', @{$blocks[$_]}) . ";T=" . join(',', @{$dc->[$_]}) . ";g=$gk->[$_][0],k=$gk->[$_][1]]" } 0..$q-1) . "\n";
              return; }
            my $sz = @{$blocks[$b]}; my $kmin = ($sz % 2) ? 1 : 0;
            for my $g (0..int(($left + 2)/2)) { for my $k ($kmin..($left + 2 - $sz - 2*$g)) { my $cont = 2*$g - 2 + $sz + $k; next if $cont < 0 || $cont > $left;
                $rec->($b + 1, $left - $cont, [@$reasons, kill_reason($blocks[$b], $g, $k, $dc->[$b], $G, $D)], [@$gk, [$g, $k]]); } } };
          $rec->(0, $n, [], []); } } }
    $badtarget++ unless $ntarget == 1;
    print "WARNING: target count $ntarget for G=$G D=$D R=(@R)\n" unless $ntarget == 1; } } }
printf "SUMMARY %s G=%d..%d: %d (G,D,R) cases; terms: %s; survivors %d; bad target counts %d; target coefficient factor -c*k = %s (non-zero needed)\n",
  $EQ, $Gmin, $Gmax, $nR, join(', ', map { "$_ $tot{$_}" } sort keys %tot), $nsurv, $badtarget, -$clin * $wt;
