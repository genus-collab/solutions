#!/usr/bin/perl
# global_suite.pl -- all variants of global_null.pl at one genus G from one set of samples (modular, p = 2^31-1).
# Variants: D in {2 (whole system: one-leaf unknowns c = 3..2G-2)} u {odd D (single-level problem of (OL))},
#           mode two (two-leaf regularity, TLR) / full (all multi-leaf T), cf 0/1 (closed-form relations), window W = G or G-1.
# Output per variant: unknowns, rank, nullity, and (single-level) whether P_{D+1} = sigma^G Xi_G(D+1) is forced to vanish.
# Usage: perl global_suite.pl G [S] [variants: e.g. "2:two:0:G 2:two:1:G-1 3:full:1:G-1 ..." or "all"]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod;
my ($G, $S, @var) = @ARGV; $S ||= 4*$G + 4; my $Lmax = 2*$G - 3;
if (!@var || $var[0] eq 'all') { @var = ();
  for my $md ('two', 'full') { for my $cf (0, 1) { push @var, "2:$md:$cf:G", "2:$md:$cf:G-1"; } }
  for (my $D = 3; $D <= $Lmax; $D += 2) { for my $md ('two', 'full') { for my $cf (0, 1) { push @var, "$D:$md:$cf:G-1"; } } } }
$| = 1;
my @samples;
for my $s (1..$S) { my $tau = md(1 + $s); my $St = kp_setup($tau, 2, $Lmax); my $rows = kp_rows($St);
  my %Pi; $Pi{$_} = $St->{Pi}->($_) for 1..2*$G; my %Ym;
  for my $T (@{$St->{allT}}) { $Ym{join(',', @$T)} = [ [@$T], $St->{Y}{join(',', @$T)} ]; }
  push @samples, { tau => $tau, sig => $St->{sig}, Pi => \%Pi, Y => \%Ym, rows => { map { ($_->[0] => $_->[1]) } @$rows } }; }
for my $v (@var) { my ($D, $mode, $cf, $Ws) = split /:/, $v; my $W = $Ws eq 'G' ? $G : $G - 1;
  my @cs = ($D+1 .. 2*$G-2); my %ci; my $nU = 0; for my $c (@cs) { for my $j (1..$W) { $ci{"$c,$j"} = $nU++; } }
  my @eqs; my %samp;
  for my $sm (@samples) { my $sig = $sm->{sig}; my $si = mi($sig); my $sinvG = 1; $sinvG = mm($sinvG, $si) for 1..$G;
    # y_L (L = 2..Lmax; index L-2 in the stored forms) as linear forms in the unknowns
    my @yform; for my $L (2..$Lmax) { my @vv = (0) x $nU; if ($L >= $D) { my $c = $L + 1; my $pw = $sinvG; for my $j (1..$W) { $pw = mm($pw, $sig); $vv[$ci{"$c,$j"}] = mm($sm->{Pi}{$c}, $pw); } } $yform[$L-2] = \@vv; }
    my $lin = sub { my $form = shift; my @vv = (0) x $nU; for my $i (0..$#$form) { next unless $form->[$i]; my $y = $yform[$i]; for my $q (0..$nU-1) { $vv[$q] = ($vv[$q] + mm($form->[$i], $y->[$q])) % $P if $y->[$q]; } } return \@vv; };
    for my $kT (keys %{$sm->{Y}}) { my ($T, $form) = @{$sm->{Y}{$kT}}; next if $mode eq 'two' && @$T != 2; my $lv = 0; $lv += $_ - 1 for @$T; next if $lv < $D;
      my $piT = 1; $piT = mm($piT, $sm->{Pi}{$_}) for @$T; my $vv = $lin->($form); my $inv = mi($piT); push @{$samp{$kT}}, [ $sig, [ map { mm($_, $inv) } @$vv ] ]; }
    if ($cf) { for my $k (sort { $a <=> $b } keys %{$sm->{rows}}) { next if $k % 2 || $k < $D; push @eqs, $lin->($sm->{rows}{$k}); } } }
  for my $kT (keys %samp) { my $L = $samp{$kT}; my $n = @$L;
    my @A = map { my $j = $_; [ map { my $sg = $L->[$_][0]; my $val = 1; if ($j >= $G) { $val = mm($val, $sg) for 1..($j-$G); } else { my $sgi = mi($sg); $val = mm($val, $sgi) for 1..($G-$j); } $val } 0..$n-1 ] } 1..$W;
    my @piv; my $row = 0;
    for my $col (0..$n-1) { my $pr; for my $rr ($row..$#A) { if ($A[$rr][$col]) { $pr = $rr; last; } } next unless defined $pr; @A[$row, $pr] = @A[$pr, $row];
      my $iv = mi($A[$row][$col]); $A[$row] = [ map { mm($iv, $_) } @{$A[$row]} ];
      for my $rr (0..$#A) { next if $rr == $row || !$A[$rr][$col]; my $f = $A[$rr][$col]; $A[$rr] = [ map { md($A[$rr][$_] - mm($f, $A[$row][$_])) } 0..$n-1 ]; }
      push @piv, $col; $row++; last if $row == @A; }
    my %isp = map { $_ => 1 } @piv;
    for my $fc (grep { !$isp{$_} } 0..$n-1) { my @w = (0) x $n; $w[$fc] = 1; for my $rr (0..$#piv) { $w[$piv[$rr]] = md(-$A[$rr][$fc]); }
      my @e = (0) x $nU; for my $s2 (0..$n-1) { next unless $w[$s2]; my $vv = $L->[$s2][1]; for my $q (0..$nU-1) { $e[$q] = ($e[$q] + mm($w[$s2], $vv->[$q])) % $P if $vv->[$q]; } }
      push @eqs, \@e; } }
  my %pivrow; my $rank = 0;
  for my $e (@eqs) { last if $rank == $nU; my @vv = @$e;
    for my $col (0..$nU-1) { next unless $vv[$col]; if (exists $pivrow{$col}) { my $f = $vv[$col]; my $pr = $pivrow{$col}; for my $q ($col..$nU-1) { $vv[$q] = md($vv[$q] - mm($f, $pr->[$q])) if $pr->[$q]; } next; }
      my $iv = mi($vv[$col]); $pivrow{$col} = [ map { mm($iv, $_) } @vv ]; $rank++; last; } }
  # forced vanishing of the lowest polynomial: rank(all) - rank(without the c = D+1 columns) == W
  my $forced = 'n/a';
  if ($rank < $nU) { my @keep = grep { $_ >= $W } 0..$nU-1; my %pr2; my $rk2 = 0;
    for my $e (@eqs) { my @vv = map { $e->[$_] } @keep; my $m = @keep;
      for my $col (0..$m-1) { next unless $vv[$col]; if (exists $pr2{$col}) { my $f = $vv[$col]; my $pr = $pr2{$col}; for my $q ($col..$m-1) { $vv[$q] = md($vv[$q] - mm($f, $pr->[$q])) if $pr->[$q]; } next; }
        my $iv = mi($vv[$col]); $pr2{$col} = [ map { mm($iv, $_) } @vv ]; $rk2++; last; } }
    $forced = ($rank - $rk2 == $W) ? 'yes' : 'no'; } else { $forced = 'yes'; }
  printf "G=%2d D=%2d mode=%-4s cf=%d W=%2d S=%d: unknowns %4d rank %4d nullity %3d  P_{D+1} forced: %s\n", $G, $D, $mode, $cf, $W, $S, $nU, $rank, $nU - $rank, $forced;
}
