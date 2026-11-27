package RefXY;
# Referee 2: Theorem 1 of the R2 report at n = 1 and n = 2, implemented from scratch (integer tau >= 1).
#   Q_m(v) = prod_{M(e2 m)} (v + 1 + k hb/e2) / prod_{M(tau m)} (v + k hb/tau) = v^m Sum_a q_a(hb) v^{-a}
#   n = 1: [x^j] F_{g,1} = j^{-2} [hb^{2g}] q_{j+1}
#   n = 2: [x^{m1} x^{m2}] F_{g,2} = (1/(m1 m2)) [hb^{2g}] [v1^-1][v2^-1] Q_{m1}(v1) Q_{m2}(v2) / ((v1-v2)^2 - K^2 hb^2),  K = (m1+m2)/2,
#          expanded in |v1| >> |v2|:  1/((v1-v2)^2 - K^2hb^2) = Sum_n (K hb)^{2n} Sum_l C(2n+1+l, l) v2^l v1^{-(2n+2+l)},
#          so  = (1/(m1 m2)) [hb^{2g}] Sum_{n,l} (K hb)^{2n} C(2n+1+l, l) q^{(m1)}_{m1-2n-1-l} q^{(m2)}_{m2+l+1}.
#   (the n = 2 sign: (-1)^{n-1} = -1 times the 2-cycle product 1/((v1-v2+K hb)(v2-v1+K hb)) = -1/((v1-v2)^2 - K^2 hb^2).)
# Dictionary P2: Xi_g(mu) = sigma^{g-1} [x^mu] F_{g,n} / prod Pi_{mu_i}.
use strict; use warnings; use Math::BigRat; use RefR2;
sub Q { Math::BigRat->new(shift) }
our %QC;
# q-coefficients: returns arrayref A[a][p] = [hb^p] q_a, a = 0..AMAX, p = 0..HB
sub qcoef { my ($tau, $m, $AMAX, $HB) = @_; my $key = "$tau|$m|$AMAX|$HB"; return $QC{$key} if exists $QC{$key};
  my $e2 = $tau + 1; my @F = map { [ map { Q(0) } 0..$HB ] } 0..$AMAX; $F[0][0] = Q(1);
  my $N2 = $e2 * $m; my $N1 = $tau * $m;
  for my $kk (1..$N2) { my $c1 = (Q($kk) - ($N2 + 1) / 2) / $e2;
    for (my $a = $AMAX; $a >= 1; $a--) { for my $p (0..$HB) { my $add = $F[$a-1][$p]->copy; $add += $c1 * $F[$a-1][$p-1] if $p >= 1; $F[$a][$p] += $add; } } }
  for my $kk (1..$N1) { my $c = (Q($kk) - ($N1 + 1) / 2) / $tau;
    for my $a (1..$AMAX) { for (my $p = $HB; $p >= 1; $p--) { $F[$a][$p] -= $c * $F[$a-1][$p-1]; } } }
  $QC{$key} = \@F; return \@F; }
sub xi1 { my ($tau, $g, $j) = @_; my $F = qcoef($tau, $j, $j + 1, 2*$g); my $sig = $tau * ($tau + 1);
  return $sig**($g-1) * $F->[$j+1][2*$g] / ($j * $j) / RefR2::Pi($j, $tau); }
sub binom { my ($n, $k) = @_; return Q(0) if $k < 0 || $k > $n; return Q(RefR2::binomI($n, $k)); }
sub F2 { my ($tau, $g, $m1, $m2) = @_; my $HB = 2*$g; my $K = Q($m1 + $m2) / 2;
  my $A = qcoef($tau, $m1, $m1 + $m2 + 2, $HB); my $B = qcoef($tau, $m2, $m1 + $m2 + 2, $HB);
  my $tot = Q(0);
  for my $n (0..$g) { my $Kp = $K**(2*$n);
    for my $l (0..$m1) { my $a = $m1 - 2*$n - 1 - $l; next if $a < 0; my $b = $m2 + $l + 1; my $c = binom(2*$n + 1 + $l, $l) * $Kp;
      # [hb^{2g}] of hb^{2n} q_a q'_b
      my $s = Q(0); for my $p (0..$HB - 2*$n) { my $x = $A->[$a][$p]; next if $x->is_zero; my $y = $B->[$b][$HB - 2*$n - $p]; next if $y->is_zero; $s += $x * $y; }
      $tot += $c * $s; } }
  return $tot / ($m1 * $m2); }
sub xi2 { my ($tau, $g, $m1, $m2) = @_; my $sig = $tau * ($tau + 1); return $sig**($g-1) * F2($tau, $g, $m1, $m2) / (RefR2::Pi($m1, $tau) * RefR2::Pi($m2, $tau)); }
1;
