#!/usr/bin/perl
# global_null.pl -- nullity of the genus-G KP system (Route A's Theorem D) with the polynomial (degree-window) ansatz,
# optionally restricted to two-leaf regularity (TLR), optionally augmented by the closed-form relations, optionally in the
# single-level setting of (OL) (one-leaf unknowns only at levels >= D).  Modular (p = 2^31-1), sampling as in Route A:
#   unknowns  P_c(sigma) = sigma^G Xi_G(c) = sum_{j=1}^{W} a_{c,j} sigma^j,  c = D+1 .. 2G-2   (W = G-1: Hurwitz window; W = G: Route A)
#   at each sample tau_s:  Y(T) = Pi_T Xi_G(T) from kpmod (KP propagation), Xi_G(T)(tau_s) = Y(T)/Pi_T as linear forms in the a's;
#   regularity:  for every multi-leaf T in the mode, (Xi_G(T)(tau_s))_s lies in span{ sigma_s^{j-G} : j = 1..W }  (left kernel);
#   closed form: relation_k(tau_s) = 0, k = D..2G-3 (the odd orders vanish identically).
# Usage: perl global_null.pl G D mode(two|full|none) cf(0|1) [W] [S]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod;
my ($G, $D, $mode, $cf, $W, $S) = @ARGV; $W //= $G - 1; $S //= 4*$G + 4; my $Lmax = 2*$G - 3;
my @cs = ($D+1 .. 2*$G-2); my %ci; my $nU = 0; for my $c (@cs) { for my $j (1..$W) { $ci{"$c,$j"} = $nU++; } }
my @eqs; my %samp;   # samp{Tkey} = [ [sigma_s, vec], ... ]
my $ninv = 0;
for my $s (1..$S) {
  my $tau = md(1 + $s); my $sig = mm($tau, md($tau + 1)); my $sinvG = mi(1); { my $si = mi($sig); $sinvG = 1; $sinvG = mm($sinvG, $si) for 1..$G; }
  my $St = kp_setup($tau, $D, $Lmax); my $Pi = $St->{Pi};
  # y_L as linear forms in the unknowns
  my @yform; for my $L ($D..$Lmax) { my @v = (0) x $nU; my $c = $L + 1; my $pw = $sinvG;
    for my $j (1..$W) { $pw = mm($pw, $sig); $v[$ci{"$c,$j"}] = mm($Pi->($c), $pw); } $yform[$L - $D] = \@v; }
  my $lin = sub { my $form = shift; my @v = (0) x $nU; for my $i (0..$#$form) { next unless $form->[$i]; my $y = $yform[$i]; for my $q (0..$nU-1) { $v[$q] = ($v[$q] + mm($form->[$i], $y->[$q])) % $P if $y->[$q]; } } return \@v; };
  if ($mode ne 'none') { for my $T (@{$St->{allT}}) { next if $mode eq 'two' && @$T != 2; my $kT = join(',', @$T);
      my $piT = 1; $piT = mm($piT, $Pi->($_)) for @$T; my $v = $lin->($St->{Y}{$kT}); my $inv = mi($piT);
      push @{$samp{$kT}}, [ $sig, [ map { mm($_, $inv) } @$v ] ]; } }
  if ($cf) { my $rows = kp_rows($St); for my $r (@$rows) { my ($k, $eq) = @$r; next if $k % 2; push @eqs, $lin->($eq); } }
}
# polynomiality: left kernel of M[s][j] = sigma_s^{j-G}, j = 1..W
for my $kT (keys %samp) { my $L = $samp{$kT}; my $n = @$L;
  my @A = map { my $j = $_; [ map { my $sg = $L->[$_][0]; my $v = mi(1); my $si = mi($sg); if ($j >= $G) { $v = 1; $v = mm($v, $sg) for 1..($j-$G); } else { $v = 1; $v = mm($v, $si) for 1..($G-$j); } $v } 0..$n-1 ] } 1..$W;
  my @piv; my $row = 0;
  for my $col (0..$n-1) { my $pr; for my $rr ($row..$#A) { if ($A[$rr][$col]) { $pr = $rr; last; } } next unless defined $pr; @A[$row, $pr] = @A[$pr, $row];
    my $iv = mi($A[$row][$col]); $A[$row] = [ map { mm($iv, $_) } @{$A[$row]} ];
    for my $rr (0..$#A) { next if $rr == $row || !$A[$rr][$col]; my $f = $A[$rr][$col]; $A[$rr] = [ map { md($A[$rr][$_] - mm($f, $A[$row][$_])) } 0..$n-1 ]; }
    push @piv, $col; $row++; last if $row == @A; }
  my %isp = map { $_ => 1 } @piv;
  for my $fc (grep { !$isp{$_} } 0..$n-1) { my @w = (0) x $n; $w[$fc] = 1; for my $rr (0..$#piv) { $w[$piv[$rr]] = md(-$A[$rr][$fc]); }
    my @e = (0) x $nU; for my $s2 (0..$n-1) { next unless $w[$s2]; my $v = $L->[$s2][1]; for my $q (0..$nU-1) { $e[$q] = ($e[$q] + mm($w[$s2], $v->[$q])) % $P if $v->[$q]; } }
    push @eqs, \@e; } }
# rank
my %pivrow; my $rank = 0;
for my $e (@eqs) { last if $rank == $nU; my @v = @$e;
  for my $col (0..$nU-1) { next unless $v[$col]; if (exists $pivrow{$col}) { my $f = $v[$col]; my $pr = $pivrow{$col}; for my $q ($col..$nU-1) { $v[$q] = md($v[$q] - mm($f, $pr->[$q])) if $pr->[$q]; } next; }
    my $iv = mi($v[$col]); $pivrow{$col} = [ map { mm($iv, $_) } @v ]; $rank++; last; } }
# which c are not fully pivoted (indicative only)
my @free; for my $c (@cs) { my $all = 1; for my $j (1..$W) { $all = 0 unless exists $pivrow{$ci{"$c,$j"}}; } push @free, $c unless $all; }
# is the lowest one-leaf polynomial forced to vanish?  test: kernel restricted to c = D+1 columns (rank with and without those columns)
printf "G=%d D=%d mode=%s cf=%d window=%d samples=%d: unknowns %d, equations %d, rank %d, nullity %d%s\n", $G, $D, $mode, $cf, $W, $S, $nU, scalar(@eqs), $rank, $nU - $rank,
  ($nU - $rank ? "; not all pivots for c = @free (indicative)" : "");
