#!/usr/bin/perl
# kp_point.pl -- the single-level problem at level D, evaluated at ONE value of the framing (exact, Math::BigRat),
# in the MV normalisation Y(T) := Pi_T Xi_G(T) (Pi_T = prod_{t in T} Pi_t), in which every relation is polynomial in tau.
# Unknowns: y_L := Y(L+1) = Pi_{L+1} Xi_G(L+1) for L = D .. Lmax (Lmax = 2G-3 for genus G).
#  KP (Theorem D, multiplied by Pi_R):  Y(R+{al,be}) = (|R|+cc) Y(R+cc) - sum_{(a,b,K) != (al,be,R)} w Pi_{R-K} E0(al-a,be-b;R-K) Y(K+{a,b}),
#      cc = al+be-1, string Y(T+1) = |T| Y(T); this defines every Y(T) as a polynomial combination of the y_L (no division).
#  Regularity at the point (first order): y_L = 0 if Pi_{L+1}(tau) = 0, Y(T) = 0 if Pi_T(tau) = 0 (mode two: |T| = 2 only; full: all T).
#  Closed-form relations, orders k = D..Lmax:  sum_T (1/|Aut T|) prod_t ((-1)^t Pi_t/t) [t^(k-D(T))] prod_t rho_t  Y(T) = 0.
# Output: rank, and whether y_D is forced to vanish at this point.
# Usage: perl kp_point.pl D Lmax tau [mode=two|full|none] [cf=1]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($D, $Lmax, $taus, $mode, $useCF) = @ARGV; $mode //= 'two'; $useCF //= 1;
my $tau = R($taus); my $sig = $tau*(1+$tau);
my $NU = $Lmax - $D + 1;
my %PiC; sub Pi { my $j = shift; return $PiC{$j} //= Pi_at($j, $tau); }
sub zero { return [ map { R(0) } 1..$NU ]; }
sub X0 { my @T = @_; my $J = 0; $J += $_ for @T; my $r = @T;
  if ($r >= 3) { my $v = R(1); $v = $v * ($J - $_) for 1..$r-3; return $v; } elsif ($r == 2) { return R(1)/$J; } else { return R(1)/($J*($J+1)); } }
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
my %Y;
for my $L ($D..$Lmax) { my $v = zero(); $v->[$L-$D] = R(1); $Y{$L+1} = $v; }
sub lvl { my $s = 0; $s += $_ - 1 for @_; return $s; }
sub YG { my @L = @_; my @Tp = sort { $b <=> $a } grep { $_ >= 2 } @L; my $nb = @L - @Tp; return zero() unless @Tp;
  my $lv = lvl(@Tp); return zero() if $lv < $D; die "YG: level $lv > Lmax" if $lv > $Lmax;
  my $J = 0; $J += $_ for @Tp; my $f = R(1); $f = $f * ($J + $_) for 0..$nb-1;
  my $v = $Y{join(',', @Tp)}; die "YG: undefined (@Tp)" unless defined $v; return [ map { $_ * $f } @$v ]; }
sub axpy { my ($y, $a, $x) = @_; return if $a->is_zero; for my $i (0..$NU-1) { $y->[$i] = $y->[$i] + $a*$x->[$i] unless $x->[$i]->is_zero; } }
my @cons; my @reg; my @cf;
my @allT;
for my $L ($D..$Lmax) {
  my @Ts = grep { @$_ >= 2 } map { [ map { $_ + 1 } @$_ ] } parts($L, $L); @Ts = sort { @$a <=> @$b } @Ts;
  my %rels;
  for my $DR (0..$L-2) { for my $q (parts($DR, $DR)) { my @Rm = map { $_ + 1 } @$q; my $w = $L - $DR + 2;
      for my $al (2..int($w/2)) { my $be = $w - $al; push @{$rels{key(@Rm, $al, $be)}}, [\@Rm, $al, $be]; } } }
  for my $T (@Ts) { my $kT = join(',', @$T); my $list = $rels{$kT} or die "no relation for $kT";
    my @sorted = sort { lvl(@{$a->[0]}) <=> lvl(@{$b->[0]}) || $b->[1] <=> $a->[1] } @$list;
    my $first = 1;
    for my $rel ($first ? @sorted : ()) { my ($Rm, $al, $be) = @$rel; my $cc = $al + $be - 1; my $JR = 0; $JR += $_ for @$Rm;
      my $val = zero(); axpy($val, R($JR + $cc), YG(@$Rm, $cc));
      for my $sb (subms($Rm)) { my ($K, $C, $w) = @$sb; my $PiC = R(1); $PiC = $PiC * Pi($_) for @$C;
        for my $a (1..$al-1) { for my $b (1..$be-1) { my $e = E0($al-$a, $be-$b, $C); next if $e->is_zero;
            my @arg = (@$K, $a, $b); my @Tp = grep { $_ >= 2 } @arg; next unless @Tp; next if lvl(@Tp) < $D;
            axpy($val, -$w * $PiC * $e, YG(@arg)); } } }
      if ($first) { $Y{$kT} = $val; $first = 0; } else { my $d = [ map { $val->[$_] - $Y{$kT}[$_] } 0..$NU-1 ]; push @cons, $d if grep { !$_->is_zero } @$d; }
      last if $ENV{FAST}; } }
  push @allT, map { [@$_] } @Ts;
}
# regularity conditions at this point
for my $L ($D..$Lmax) { push @reg, $Y{$L+1} if Pi($L+1)->is_zero; }
if ($mode ne 'none') { for my $T (@allT) { next if $mode eq 'two' && @$T != 2; my $p = R(1); $p = $p * Pi($_) for @$T; push @reg, $Y{join(',', @$T)} if $p->is_zero; } }
# closed-form relations
if ($useCF) {
  my $M = $Lmax - $D + 1;
  my $Tt = T_of_t_series($M+1); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $M+1), $Tt, $M+1);
  my %rho;
  for my $j (2..$Lmax+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $M+1), $Tt, $M+1);
    my $p = [ R(1), map { R(0) } 1..$M+1 ]; $p = series_mul($p, $h1, $M+1) for 1..$j; $rho{$j} = series_mul($hj, series_inv($p, $M+1), $M+1); }
  for my $k ($D..$Lmax) { my $eq = zero();
    for my $LT ($D..$k) { for my $q (parts($LT, $LT)) { my @T = map { $_ + 1 } @$q; my $top = R(1); my $pr = [ R(1), map { R(0) } 1..$M ];
        for my $t (@T) { $top = $top * R((-1)**$t) * Pi($t) / $t; $pr = series_mul($pr, $rho{$t}, $M); }
        my $co = $top * $pr->[$k - $LT] / aut(@T); next if $co->is_zero; axpy($eq, $co, $Y{join(',', @T)}); } }
    push @cf, [$k, $eq]; } }
sub rowspace { my @rows = map { [ @$_ ] } @_; my $rk = 0;
  for my $col (0..$NU-1) { my $pr; for my $i ($rk..$#rows) { unless ($rows[$i][$col]->is_zero) { $pr = $i; last; } } next unless defined $pr;
    @rows[$rk, $pr] = @rows[$pr, $rk]; my $inv = R(1)/$rows[$rk][$col]; $rows[$rk] = [ map { $_ * $inv } @{$rows[$rk]} ];
    for my $i (0..$#rows) { next if $i == $rk || $rows[$i][$col]->is_zero; my $f = $rows[$i][$col]; for my $j ($col..$NU-1) { $rows[$i][$j] = $rows[$i][$j] - $f*$rows[$rk][$j]; } }
    $rk++; }
  return ($rk, [ @rows[0..$rk-1] ]); }
sub forced { my ($rows, $col) = @_; for my $r (@$rows) { my $nz = grep { !$_->is_zero } @$r; return 1 if $nz == 1 && !$r->[$col]->is_zero; } return 0; }
my $cfodd = 0; for my $c (@cf) { $cfodd++ if ($c->[0] % 2 == 1) && grep { !$_->is_zero } @{$c->[1]}; }
if ($ENV{DUMPY}) { for my $k (sort keys %Y) { print "Y $k : ", join(" ", map { "$_" } @{$Y{$k}}), "
"; } exit 0; }
my ($r0) = rowspace(@cons);
my ($r1, $rs1) = rowspace(@reg);
my ($r2, $rs2) = rowspace(@reg, map { $_->[1] } @cf);
my ($r3) = rowspace(map { $_->[1] } @cf);
if ($ENV{PRINTCF}) { for my $c (@cf) { print "REL $c->[0]: ", join(" ", map { "$_" } @{$c->[1]}), "
"; } }
printf "D=%d Lmax=%d tau=%s mode=%s: %d unknowns y_%d..y_%d; KP consistency rank %d; regularity eqs %d (rank %d, y_D forced %s); closed form rank %d (odd orders non-zero: %d); together rank %d, y_D forced: %s\n",
  $D, $Lmax, $taus, $mode, $NU, $D, $Lmax, $r0, scalar(@reg), $r1, forced($rs1,0) ? 'yes' : 'no', $r3, $cfodd, $r2, forced($rs2,0) ? 'YES' : 'no';
