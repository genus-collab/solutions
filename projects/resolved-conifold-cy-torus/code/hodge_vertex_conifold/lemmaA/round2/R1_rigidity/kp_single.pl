#!/usr/bin/perl
# kp_single.pl -- the inductive single-level problem at level D, pointwise in the framing tau (exact, Math::BigRat).
# Hypotheses (those of (OL)): Lemma A in genera < G, Xi_G(T) = 0 for D(T) < D, and the levels used are in the range (<= 2G-3).
# Unknowns: one-leaf values v_c = Xi_G(c), c = D+1 .. Lmax+1 (numbers at the fixed tau).
# (i)  KP (Route A, Theorem D), every relation (R, al, be), al <= be, level L = D(R)+al+be-2 in [D, Lmax]:
#        sum_{a<=al, b<=be} sum_{K subset R} Pi_a Pi_b Xi(K+{a,b}) E0(al-a, be-b; R-K) = Pi_{al+be-1} Xi(R+{1, al+be-1}),
#      string for unit parts, E0 = coefficients of exp(sum z^-a w^-b sigma Pi_a Pi_b Xi_0(.+{a,b})) (labelled), Xi_0(T) = Gamma(J)/Gamma(J-r+3).
#      One relation per multi-part T defines Xi(T) (two largest parts merged); all other relations are recorded as consistency equations.
# (ii) closed-form relations (single-vertex graph at polar order k, brief s.2(e)): for k = D .. Lmax
#        [t^k] sum_T (1/|Aut T|) prod_{t in T} b_t(t) Xi(T) = 0,   b_j = (-1)^j (Pi_j^2/j) t^(j-1) H_j(T)/H_1(T)^j,  H_j = 2F1(-j tau,-j tau; j; 1-e^-t).
# Output: rank of the consistency equations alone, and of everything; whether v_{D+1} is forced to vanish.
# Usage: perl kp_single.pl D Lmax tau [closed=1]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($D, $Lmax, $taus, $useCF) = @ARGV; $useCF //= 1;
my $tau = R($taus); my $sig = $tau*(1+$tau);
my $NU = $Lmax - $D + 1;            # unknowns v_{D+1} .. v_{Lmax+1}
my %PiC; sub Pi { my $j = shift; return $PiC{$j} //= Pi_at($j, $tau); }
sub zero { return [ map { R(0) } 1..$NU ]; }
sub X0 { my @T = @_; my $J = 0; $J += $_ for @T; my $r = @T;
  if ($r >= 3) { my $v = R(1); $v = $v * ($J - $_) for 1..$r-3; return $v; } elsif ($r == 2) { return R(1)/$J; } else { return R(1)/($J*($J+1)); } }
# sub-multisets of a multiset (desc array) with labelled weights: returns list of [K, C, w]
my %SubC;
sub subms { my $S = shift; my $k = join(',', @$S); return @{ $SubC{$k} //= do {
  my %m; $m{$_}++ for @$S; my @v = sort { $b <=> $a } keys %m; my @out = ([[], [], R(1)]);
  for my $val (@v) { my $n = $m{$val}; my @new;
    for my $o (@out) { for my $i (0..$n) { push @new, [ [@{$o->[0]}, ($val) x $i], [@{$o->[1]}, ($val) x ($n-$i)], $o->[2]*binom($n,$i) ]; } }
    @out = @new; }
  \@out } }; }
my %E0C;
sub E0 { my ($m, $n, $S) = @_; my $k = "$m;$n;" . join(',', @$S); return $E0C{$k} if exists $E0C{$k};
  my $res;
  if ($m == 0 && $n == 0) { $res = @$S ? R(0) : R(1); }
  elsif ($m <= 0 || $n <= 0) { $res = R(0); }
  else { $res = R(0);
    for my $sb (subms($S)) { my ($K, $C, $w) = @$sb;
      for my $a1 (1..$m) { for my $b1 (1..$n) { my $e = E0($m-$a1, $n-$b1, $C); next if $e->is_zero;
          $res = $res + $w * $a1 * $sig * Pi($a1) * Pi($b1) * X0(@$K, $a1, $b1) * $e; } } }
    $res = $res / $m; }
  $E0C{$k} = $res; return $res; }
# genus-G values as linear forms
my %XI;
for my $c ($D+1..$Lmax+1) { my $v = zero(); $v->[$c-$D-1] = R(1); $XI{$c} = $v; }
sub lvl { my $s = 0; $s += $_ - 1 for @_; return $s; }
sub XG { my @L = @_; my @Tp = sort { $b <=> $a } grep { $_ >= 2 } @L; my $nb = @L - @Tp; return zero() unless @Tp;
  my $lv = lvl(@Tp); return zero() if $lv < $D; die "XG: level $lv > Lmax" if $lv > $Lmax;
  my $J = 0; $J += $_ for @Tp; my $f = R(1); $f = $f * ($J + $_) for 0..$nb-1;
  my $v = $XI{join(',', @Tp)}; die "XG: undefined value for (@Tp)" unless defined $v; return [ map { $_ * $f } @$v ]; }
sub axpy { my ($y, $a, $x) = @_; for my $i (0..$NU-1) { $y->[$i] = $y->[$i] + $a*$x->[$i] unless $x->[$i]->is_zero; } }
my @eqs; my @cfeqs;
for my $L ($D..$Lmax) {
  # all T at level L with r >= 2 parts, fewer parts first
  my @Ts = grep { @$_ >= 2 } map { [ map { $_ + 1 } @$_ ] } parts($L, $L); @Ts = sort { @$a <=> @$b } @Ts;
  # relations with main term at this level, grouped by main T
  my %rels;
  for my $DR (0..$L-2) { for my $q (parts($DR, $DR)) { my @Rm = map { $_ + 1 } @$q; my $w = $L - $DR + 2;   # al + be
      for my $al (2..int($w/2)) { my $be = $w - $al; push @{$rels{key(@Rm, $al, $be)}}, [\@Rm, $al, $be]; } } }
  for my $T (@Ts) { my $kT = join(',', @$T); my $list = $rels{$kT} or die "no relation for $kT";
    # canonical first: largest al+be (smallest D(R)), then larger al
    my @sorted = sort { lvl(@{$a->[0]}) <=> lvl(@{$b->[0]}) || $b->[1] <=> $a->[1] } @$list;
    my $first = 1;
    for my $rel (@sorted) { my ($Rm, $al, $be) = @$rel; my $cc = $al + $be - 1;
      my $rest = zero();
      for my $sb (subms($Rm)) { my ($K, $C, $w) = @$sb;
        for my $a (1..$al-1) { for my $b (1..$be-1) { my $e = E0($al-$a, $be-$b, $C); next if $e->is_zero;
            my @arg = (@$K, $a, $b); my @Tp = grep { $_ >= 2 } @arg; next unless @Tp; next if lvl(@Tp) < $D;
            axpy($rest, $w * Pi($a) * Pi($b) * $e, XG(@arg)); } } }
      axpy($rest, -Pi($cc), XG(@$Rm, 1, $cc));
      my $mc = Pi($al) * Pi($be);
      if ($first) { $XI{$kT} = [ map { -$_ / $mc } @$rest ]; $first = 0; }
      else { my $eq = [ map { $mc * $XI{$kT}[$_] + $rest->[$_] } 0..$NU-1 ]; push @eqs, $eq; } } } }
# closed-form relations
if ($useCF) {
  my $M = $Lmax - $D + 1;   # relative orders needed
  my $Tt = T_of_t_series($M+1); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $M+1), $Tt, $M+1);
  my %rho;
  for my $j (2..$Lmax+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $M+1), $Tt, $M+1);
    my $p = [ R(1), map { R(0) } 1..$M+1 ]; $p = series_mul($p, $h1, $M+1) for 1..$j; $rho{$j} = series_mul($hj, series_inv($p, $M+1), $M+1); }
  for my $k ($D..$Lmax) { my $eq = zero();
    for my $LT ($D..$k) { for my $q (parts($LT, $LT)) { my @T = map { $_ + 1 } @$q; my $top = R(1); my $pr = [ R(1), map { R(0) } 1..$M ];
        for my $t (@T) { $top = $top * R((-1)**$t) * Pi($t) * Pi($t) / $t; $pr = series_mul($pr, $rho{$t}, $M); }
        my $co = $top * $pr->[$k - $LT] / aut(@T); next if $co->is_zero; axpy($eq, $co, $XI{join(',', @T)}); } }
    push @cfeqs, [$k, $eq]; } }
# ranks
sub rank_of { my @rows = map { [ @$_ ] } @_; my $rk = 0; my @piv;
  for my $col (0..$NU-1) { my $pr; for my $i ($rk..$#rows) { unless ($rows[$i][$col]->is_zero) { $pr = $i; last; } } next unless defined $pr;
    @rows[$rk, $pr] = @rows[$pr, $rk]; my $inv = R(1)/$rows[$rk][$col];
    for my $i (0..$#rows) { next if $i == $rk || $rows[$i][$col]->is_zero; my $f = $rows[$i][$col] * $inv; for my $j ($col..$NU-1) { $rows[$i][$j] = $rows[$i][$j] - $f*$rows[$rk][$j]; } }
    push @piv, $col; $rk++; }
  # is v_{D+1} (col 0) forced to vanish? i.e. is e_0 in the row space
  my @rs = @rows[0..$rk-1]; my $forced = 0;
  for my $r (@rs) { my $nz = 0; for my $j (0..$NU-1) { $nz++ unless $r->[$j]->is_zero; } $forced = 1 if $nz == 1 && !$r->[0]->is_zero; }
  return ($rk, $forced); }
my ($rk1, $f1) = rank_of(@eqs);
my $nonzero1 = grep { my $r = $_; grep { !$_->is_zero } @$r } @eqs;
printf "D=%d Lmax=%d tau=%s: unknowns v_%d..v_%d (%d); KP consistency equations %d (non-zero %d), rank %d, v_%d forced: %s\n",
  $D, $Lmax, $taus, $D+1, $Lmax+1, $NU, scalar(@eqs), $nonzero1, $rk1, $D+1, $f1 ? "yes" : "no";
if ($useCF) {
  for my $ce (@cfeqs) { my ($k, $eq) = @$ce; print "  closed-form order $k: ", join("  ", map { "v_" . ($D+1+$_) . ":" . $eq->[$_] } grep { !$eq->[$_]->is_zero } 0..$NU-1), "\n"; }
  my ($rk2, $f2) = rank_of(@eqs, map { $_->[1] } @cfeqs);
  printf "  with closed-form relations (orders %d..%d): rank %d of %d, v_%d forced: %s\n", $D, $Lmax, $rk2, $NU, $D+1, $f2 ? "yes" : "no";
}
