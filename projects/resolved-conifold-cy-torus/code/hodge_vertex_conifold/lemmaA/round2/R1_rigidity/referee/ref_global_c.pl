#!/usr/bin/perl
# ref_global.pl -- referee's independent recomputation of the "global nullities" (route s.12), whole system (D = 2):
# unknowns a_{c,k}: sigma^G Xi_G(c) = Sum_{k=1}^{W} a_{c,k} sigma^k, c = 3..2G-2 (Xi_G(2) = 0 by Theorem A).
# Constraints, at S sample framings tau_s = s+1 (mod p):
#   regularity: for every multi-leaf T in the mode (two: |T| = 2, a+b <= 2G-1; full: every T at level <= 2G-3),
#     Xi_G(T)(tau_s) = Y(T)(tau_s)/Pi_T(tau_s) must lie in span_s{ sigma_s^(k-G) : k = 1..W } (eliminated by left kernels);
#   cf: relation_k = 0 for even k <= 2G-3 (input I1; odd k are automatic).
# Y(T) from refkp's propagation.  Prints unknowns, rank, nullity.
# Usage: perl ref_global.pl G mode(two|full) cf(0|1) W(G|G-1) [S]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use refkp;
my ($G, $mode, $cf, $Ws, $S) = @ARGV; set_mode('mod'); set_prime($ENV{PRIME}) if $ENV{PRIME};
my $W = $Ws eq 'G' ? $G : $G - 1; $S ||= 4*$G + 6; my $K = 2*$G - 3;
my @cs = (($ENV{CMIN} // 3)..2*$G-2); my %ci; my $nU = 0; for my $c (@cs) { for my $k (1..$W) { $ci{"$c,$k"} = $nU++; } }
my %samp; my @eqs;
for my $s (1..$S) { my $tau = NN($s + 1); my $St = kp_propagate($tau, $K); my $Pi = $St->{Pi}; my $sig = $St->{sig};
  my $sgi = IV($sig); my $sG = ONE(); $sG = ML($sG, $sgi) for 1..$G;
  # y_L as linear forms in a
  my @yf; for my $L (1..$K) { my $v = vzero($nU); my $c = $L + 1; if ($c >= 3) { my $pw = $sG; for my $k (1..$W) { $pw = ML($pw, $sig); $v->[$ci{"$c,$k"}] = ML($Pi->($c), $pw); } } $yf[$L] = $v; }
  my $lin = sub { my $form = shift; my $v = vzero($nU); for my $L (1..$K) { vaxpy($v, $form->[$L-1], $yf[$L]); } return $v; };
  for my $Tk (keys %{$St->{Y}}) { my @T = split /,/, $Tk; next if @T < 2; next if $mode eq 'two' && @T != 2;
    my $piT = ONE(); $piT = ML($piT, $Pi->($_)) for @T; my $v = vscale(IV($piT), $lin->($St->{Y}{$Tk})); push @{$samp{$Tk}}, [ $sig, $v ]; }
  if ($cf) { my $Tt = Tser($K); my $h1 = ser_comp(hyp2f1(NG($tau), NG($tau), ONE(), $K), $Tt, $K); my %rho;
    for my $j (2..$K+1) { my $hj = ser_comp(hyp2f1(NG(ML(NN($j), $tau)), NG(ML(NN($j), $tau)), NN($j), $K), $Tt, $K);
      my $pw = [ ONE(), map { ZERO() } 1..$K ]; $pw = ser_mul($pw, $h1, $K) for 1..$j; $rho{$j} = ser_mul($hj, ser_inv($pw, $K), $K); }
    my @rel = map { vzero($K) } 0..$K;
    for my $Dl (1..$K) { for my $T (parts_ge2($Dl, $Dl + 1)) { my $top = ONE(); my $pr = [ ONE(), map { ZERO() } 1..$K ];
        for my $t (@$T) { $top = ML($top, ML(NN(($t % 2) ? -1 : 1), ML($Pi->($t), IV(NN($t))))); $pr = ser_mul($pr, $rho{$t}, $K - $Dl); }
        $top = ML($top, IV(aut_f(@$T))); my $v = $St->{Y}{mkey(@$T)}; for my $k ($Dl..$K) { vaxpy($rel[$k], ML($top, $pr->[$k - $Dl]), $v); } } }
    for (my $k = 2; $k <= $K; $k += 2) { push @eqs, $lin->($rel[$k]); } } }
# eliminate the free polynomials of the multi-leaf values: left kernel of M_{s,k} = sigma_s^(k-G)
for my $Tk (keys %samp) { my $L = $samp{$Tk}; my $n = @$L;
  my @A = map { my $k = $_; [ map { my $sg = $L->[$_][0]; my $v = ONE(); if ($k >= $G) { $v = ML($v, $sg) for 1..($k-$G); } else { my $si = IV($sg); $v = ML($v, $si) for 1..($G-$k); } $v } 0..$n-1 ] } 1..$W;
  # row-reduce A (W x n); kernel vectors w with A w = 0 give equations Sum_s w_s X_s = 0
  my @piv; my $row = 0;
  for my $col (0..$n-1) { my $pr; for my $r ($row..$#A) { unless (IZ($A[$r][$col])) { $pr = $r; last; } } next unless defined $pr; @A[$row, $pr] = @A[$pr, $row];
    my $iv = IV($A[$row][$col]); $A[$row] = [ map { ML($iv, $_) } @{$A[$row]} ];
    for my $r (0..$#A) { next if $r == $row || IZ($A[$r][$col]); my $f = $A[$r][$col]; $A[$r] = [ map { SB($A[$r][$_], ML($f, $A[$row][$_])) } 0..$n-1 ]; }
    push @piv, $col; $row++; last if $row == @A; }
  my %isp = map { $_ => 1 } @piv;
  for my $fc (grep { !$isp{$_} } 0..$n-1) { my @w = map { ZERO() } 1..$n; $w[$fc] = ONE(); for my $r (0..$#piv) { $w[$piv[$r]] = NG($A[$r][$fc]); }
    my $e = vzero($nU); for my $s2 (0..$n-1) { vaxpy($e, $w[$s2], $L->[$s2][1]); } push @eqs, $e; } }
# rank
my %pr; my $rank = 0;
for my $e (@eqs) { last if $rank == $nU; my @v = @$e;
  for my $col (0..$nU-1) { next if IZ($v[$col]); if (exists $pr{$col}) { my $f = $v[$col]; my $p = $pr{$col}; for my $q ($col..$nU-1) { $v[$q] = SB($v[$q], ML($f, $p->[$q])); } next; }
    my $iv = IV($v[$col]); $pr{$col} = [ map { ML($iv, $_) } @v ]; $rank++; last; } }
printf "ref_global G=%d mode=%s cf=%d W=%d S=%d: unknowns %d, equations %d, rank %d, nullity %d\n", $G, $mode, $cf, $W, $S, $nU, scalar(@eqs), $rank, $nU - $rank;
