# ref_claimE.pl : referee checks of the proof of Claim E (report s.2), SYMBOLIC where the author checked numerically.
#  (V1) X(w/(1+w)) - X(w) = -log(1+w) + w/2 + w/(2(1+w)), X(w) = Sum_m chi_m w^{2m}, chi_m = -B_{2m}/(2m)   (to order w^NW)
#  (V2) Sum_{M>=1} C(M;a,0) w^{2M+2} = X(w)[(a-1)a w^2 + 2 - (1+w)^a - (1-w)^a] + (2a-1)w^2
#                                        - (w/2)[(1+w)^a - (1-w)^a + (1+w)^{a-1} - (1-w)^{a-1}]
#       as an identity of power series whose coefficients are POLYNOMIALS IN a (the n-dependence is the trivial term 2n(a-1) w^2 X),
#       C(M;a,0) = chi_M 2M(2a-2M-1) - 2 d/da binom(a,2M+2) - 2 Sum_{m=1}^{M-1} binom(a-2m, 2M+2-2m) chi_m.
#  (P)  P_j(w) = 0 as an identity of polynomials in j, coefficient of w^{2N}, N <= NP.
#  (R)  R_j(w) has degree <= j+1 for integers 1 <= j <= JR.
#  (E)  E(h,j) = 0 directly from its definition, 1 <= j <= 2h-1, h <= HE.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RL;
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($NW, $NP, $JR, $HE) = @ARGV; $NW //= 32; $NP //= 18; $JR //= 30; $HE //= 20;
my %chi; $chi{$_} = -RL::bern(2*$_) / (2*$_) for 1..int($NW/2) + 2;
# ---------- (V1) ----------
{ my $W = $NW; my @X = map { Q(0) } 0..$W; $X[2*$_] = $chi{$_} for 1..int($W/2);
  my @u = map { Q(0) } 0..$W; $u[$_] = Q((-1) ** ($_ - 1)) for 1..$W;             # w/(1+w)
  my @comp = map { Q(0) } 0..$W; my @p = map { Q(0) } 0..$W; $p[0] = Q(1);
  for my $k (1..$W) { my @np = map { Q(0) } 0..$W; for my $i (0..$W) { next if $p[$i] == 0; for my $l (1..$W-$i) { $np[$i+$l] += $p[$i] * $u[$l] } } @p = @np; next unless $k % 2 == 0; for my $i (0..$W) { $comp[$i] += $X[$k] * $p[$i] } }
  my $bad = 0; for my $i (0..$W) { my $lhs = $comp[$i] - $X[$i]; my $rhs = Q(0); if ($i >= 1) { $rhs += -Q((-1) ** ($i - 1)) / $i + Q((-1) ** ($i - 1)) / 2 } $rhs += Q(1)/2 if $i == 1; $bad++ if $lhs != $rhs }
  print "(V1) to order w^$W: ", ($bad ? "$bad FAILURES" : "holds"), "\n"; }
# ---------- (V2) symbolic in a ----------
{ my $W = $NW; my $binA = sub { my ($c, $k) = @_; return RL::psc(RL::ffj($c, $k), Q(1) / RL::fact($k)) };   # binom(a + c, k) as polynomial in a
  my @lhs = map { [Q(0)] } 0..$W;
  for my $M (1..int(($W - 2)/2)) { my $C = RL::psc([Q(-2*$M - 1) * 2*$M, Q(4*$M)], $chi{$M});               # chi_M 2M(2a-2M-1)
    $C = RL::padd($C, RL::pder($binA->(0, 2*$M + 2)), Q(-2));
    for my $m (1..$M-1) { $C = RL::padd($C, RL::psc($binA->(-2*$m, 2*$M + 2 - 2*$m), $chi{$m}), Q(-2)) }
    $lhs[2*$M + 2] = $C; }
  my @rhs = map { [Q(0)] } 0..$W;
  # bracket b(w) = (a-1)a w^2 + 2 - (1+w)^a - (1-w)^a
  my @br = map { [Q(0)] } 0..$W; $br[0] = [Q(2)]; $br[2] = [Q(0), Q(-1), Q(1)];
  for my $k (0..$W) { my $bk = $binA->(0, $k); $br[$k] = RL::padd($br[$k], $bk, Q(-1)); $br[$k] = RL::padd($br[$k], $bk, Q(-1) * Q((-1) ** $k)) }
  for my $i (0..$W) { for my $m (1..int($i/2)) { $rhs[$i] = RL::padd($rhs[$i], RL::psc($br[$i - 2*$m], $chi{$m})) } }   # X(w) * b(w)
  $rhs[2] = RL::padd($rhs[2], [Q(-1), Q(2)]);
  for my $k (0..$W-1) { my $d = RL::padd(RL::padd(RL::padd($binA->(0, $k), $binA->(0, $k), Q(-1) * Q((-1) ** $k)), $binA->(-1, $k)), $binA->(-1, $k), Q(-1) * Q((-1) ** $k));
    $rhs[$k + 1] = RL::padd($rhs[$k + 1], $d, Q(-1)/2) }
  my $bad = 0; for my $i (0..$W) { $bad++ unless RL::peq($lhs[$i], $rhs[$i]) }
  print "(V2) as polynomials in a, to order w^$W: ", ($bad ? "$bad FAILURES" : "holds"), "\n"; }
# ---------- (P) symbolic in j ----------
$RL::NT = 2*$NP + 4; my $L = RL::logH(); my $one = RL::szero(); $one->[0] = Q(1);
my %T; for my $n (0..$NP) { $T{$n} = RL::pmul(RL::ffj(-1, 2*$n), RL::coeff_expj($one, $L, 2*$n)) }
{ my $bad = 0;
  for my $N (1..$NP) { my $c = RL::psc(RL::pmul($T{$N - 1}, [Q(-2*$N + 1), Q(1)]), Q(2*($N - 1)));
    for my $n (0..$N-2) { $c = RL::padd($c, RL::pmul($T{$n}, RL::psc(RL::ffj(-2*$n, 2*$N - 2*$n), Q(1) / RL::fact(2*$N - 2*$n))), Q(-2)) }
    $bad++ unless RL::pzero($c) }
  print "(P) P_j(w) = 0: coefficients of w^{2N}, N = 1..$NP, vanish as polynomials in j: ", ($bad ? "$bad FAILURES" : "yes"), "\n"; }
# ---------- (R) degree bound for integer j ----------
{ my $bad = 0; my $W = 2*$JR + 6; my %Tv;
  for my $j (1..$JR) { my @R = map { Q(0) } 0..$W;
    for (my $n = 0; 2*$n < $j; $n++) { my $a = $j - 2*$n; my $tn = RL::pev($T{$n} // do { $T{$n} = RL::pmul(RL::ffj(-1, 2*$n), RL::coeff_expj($one, $L, 2*$n)) }, Q($j));
      my @br = map { Q(0) } 0..$W; $br[2] += 2*$a - 1;
      for my $k (0..$W-1) { my $d = RL::binq($a, $k) * (1 - (-1) ** $k) + RL::binq($a - 1, $k) * (1 - (-1) ** $k); $br[$k + 1] -= $d / 2 }
      for my $i (0..$W - 2*$n) { $R[$i + 2*$n] += $tn * $br[$i] } }
    my @nz = grep { $R[$_] != 0 } 0..$W; my $deg = @nz ? $nz[-1] : -1; $bad++ if $deg > $j + 1; }
  print "(R) deg R_j <= j+1 for j = 1..$JR: ", ($bad ? "$bad FAILURES" : "yes"), "\n"; }
# ---------- (E) direct ----------
{ $RL::NT = 2*$HE + 4; my $L2 = RL::logH(); my $one2 = RL::szero(); $one2->[0] = Q(1);
  my %Tp; for my $n (0..$HE) { $Tp{$n} = RL::pmul(RL::ffj(-1, 2*$n), RL::coeff_expj($one2, $L2, 2*$n)) }
  my %ch; $ch{$_} = -RL::bern(2*$_) / (2*$_) for 1..$HE + 1; my %ze; $ze{$_} = (2*$_ - 1) * RL::bern(2*$_) / (2*$_) for 1..$HE + 1;
  my $c = sub { Q(2) / RL::fact(2 * $_[0]) };
  my $ff = sub { my ($x, $n) = @_; my $p = Q(1); $p *= ($x - $_) for 0..$n-1; $p };
  my $dff = sub { my ($x, $n) = @_; RL::pev(RL::pder(RL::ffj(0, $n)), $x) };
  my ($tot, $bad) = (0, 0);
  for my $h (1..$HE) { for my $j (1..2*$h - 1) { my $jq = Q($j); my $E = Q(0);
      for (my $n = 0; 2*$n < $j; $n++) { my $tn = RL::pev($Tp{$n}, $jq);
        my $br = (2*$h + 1) * ($j - 2*$h) * $ch{$h - $n} - $c->($h + 1 - $n) * $dff->($jq - 2*$n, 2*$h + 2 - 2*$n) - $jq * $ze{$h - $n};
        for my $hp (1..$h-1) { next unless 2*$hp >= $j + 1; $br -= $c->($h + 1 - $hp) * $ff->($jq - 2*$hp, 2*$h + 2 - 2*$hp) * $ch{$hp - $n} }
        $E += $tn * $br; }
      $tot++; $bad++ if $E != 0; } }
  print "(E) E(h,j) = 0 for 1 <= j <= 2h-1, h <= $HE: $tot pairs, ", ($bad ? "$bad FAILURES" : "all zero"), "\n"; }
