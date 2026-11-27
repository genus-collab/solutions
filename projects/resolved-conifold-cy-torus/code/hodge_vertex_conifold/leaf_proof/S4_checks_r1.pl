#!/usr/bin/perl
# S4_checks_r1.pl -- exact spot checks (Perl + Math::BigRat; no Mathematica) for the parts of leaf_proof/S4_legs_closed.md
# that were added or rewritten in revision 1 (referee report of 2026-10-31). Nothing here uses a closed form of the rows:
# the rows at U_2 are built from the code's definitions: the Wick row wickZ of wickDirect in W_k = (U d/dU)^k spotU at U_2(t),
# computed as the Gaussian moments (D1) of S1 s.2.2; the frame-B twist h_w^-1 of twistRows[hCoeffs[{c,-x,-1},-1]]; the H-row of
# rmatrixCore; Delta_i = 1/etacan[[i,i]] as the literal quadratic form of rc_lib.m s.1; EdN of evalQ.
#  Part A (fixed t != 0, exact rational points (x, tau) with rational delta):
#   A1  roots: Q x U^2 + Q (1-x) U - 1 = 0 at U_2 = 1 - tau s_2, s_2 = (c + delta)/(2 x tau); W_1(U_2) = 0;
#   A2  Step A of Lemma 10.5: the code's Wick row at U_2 (coordinate y = log U, amplitude 1) equals
#       U_2 Lnorm[phi_tau, (1 - tau s)^-1; s_2] for z^0..z^5;
#   A3  W_2(U_2) = U_2^2 phi''_tau(s_2)/tau^2 = x delta s*^2,  phi''_tau(s_2) = delta/(x s_2^2 U_2^2),
#       Delta_2 (quadratic form) = x^2 W_2(U_2) = x^3 delta s*^2   (V1 at U_2, tensor rule);
#   A4  (1 - tau) s_2/(1 - tau s_2) = 1 - x s_2 = L_2 - x - c (1 - tau)/tau   (d f0/dt = L_2 - mu, Step B).
#  Part B (t-series at rational x, coefficients t^0..t^MAXC):
#   B1  roots U_1, U_2 of LeqU; W_1(U_2) = 0; L_2 = x s*; Delta_2 (quadratic form) = x^3 delta s*^2 = x^2 W_2(U_2)
#       = x^2 U_2^2 phi''_tau(s_2)/tau^2; Delta_2'(0) = 0;
#   B2  d/dt [phi_tau(s_2) + c log tau] = L_2 - x   (Step B of Lemma 10.5);
#   B3  Row 2 at U_2 (Lemma 10.5): t (D_2 r_{1,2})[m-1] = (t L_1)(r_{1,2}[m] - L_2 r_{0,2}[m]), m = 1..MZ;
#   B4  Prop. 10.6(a): r^pt_{0,2}(z; 0) = h_w^-1, r_{0,2}(z; 0) = h_w^-2;  10.6(b): r_{1,2}(z; 0) = (x/c) r_{0,2}(z; 0);
#       [t^1] r_{0,2} = 0 (Remark 10.6');
#   B5  E_22 at t = 0 from the code's edge formula (1.1) against 10.6(c), (x^3/c) [1 - h_w(z)^2 h_w(w)^2]/(z + w);
#       E_22(0,0) = -x^2 (1 + x + x^2)/(6 c^2); the 2x2 minor (x^6/c^2) S_1 (72 S_2 + 5 S_1^3)/77760 != 0; and the number of
#       entries at which the value -(x^3/2)/c^(m+n+2) printed in conifold_regularity_note.md F3 agrees;
#   B6  CF.2(d): rho = (t^2 Delta_1/(-x c^3))^(-1/2) (quadratic form, branch rho(0) = 1) equals
#       (tau/t) 2/((1 - u)^(1/4) (1 + (1 - u)^(1/2))), u = 4 x tau/c^2 (no branch ambiguity, also for c < 0);
#   B7  Remark 10.6'': <v_2(z), v_2(-z)> = Delta_2 (K_2 = 1), z^0..z^MZ.
# Usage: perl S4_checks_r1.pl [x1 x2 ...]      (default: 3/7 5/3 -5/2 2/9)
use strict; use warnings; use Math::BigRat;

sub rat { my $v = shift; return ref($v) ? $v->copy : Math::BigRat->new("$v"); }
our $NT = 9;           # t-series kept to t^NT
my $MAXC = 6;          # t-orders compared (derivatives and divisions by t lose at most 3 orders)
my $MZ = 6;            # z-orders of the rows
my @XS = @ARGV ? @ARGV : ("3/7", "5/3", "-5/2", "2/9");

sub szero { return [ map { Math::BigRat->new(0) } 0..$NT ]; }
sub cst { my $a = shift; my $s = szero(); $s->[0] = rat($a); return $s; }
sub sadd { my ($a, $b) = @_; return [ map { $a->[$_] + $b->[$_] } 0..$NT ]; }
sub ssub { my ($a, $b) = @_; return [ map { $a->[$_] - $b->[$_] } 0..$NT ]; }
sub ssc  { my ($k, $a) = @_; my $kk = rat($k); return [ map { $a->[$_] * $kk } 0..$NT ]; }
sub smul { my ($a, $b) = @_; my $r = szero();
  for my $i (0..$NT) { next if $a->[$i]->is_zero; for my $j (0..$NT-$i) { next if $b->[$j]->is_zero; $r->[$i+$j] += $a->[$i] * $b->[$j]; } }
  return $r; }
sub sinv { my $a = shift; die "sinv: zero constant term" if $a->[0]->is_zero; my $r = szero(); $r->[0] = rat(1) / $a->[0];
  for my $n (1..$NT) { my $s = rat(0); for my $k (1..$n) { $s += $a->[$k] * $r->[$n-$k]; } $r->[$n] = (rat(0) - $s) / $a->[0]; }
  return $r; }
sub sdiv { my ($a, $b) = @_; return smul($a, sinv($b)); }
sub sder { my $a = shift; my $r = szero(); for my $n (1..$NT) { $r->[$n-1] = $a->[$n] * $n; } return $r; }      # loses t^NT
sub sint { my $a = shift; my $r = szero(); for my $n (0..$NT-1) { $r->[$n+1] = $a->[$n] / ($n + 1); } return $r; }
sub sdownt { my $a = shift; die "sdownt: nonzero constant term" unless $a->[0]->is_zero;                          # a/t, loses t^NT
  my $r = szero(); for my $n (0..$NT-1) { $r->[$n] = $a->[$n+1]->copy; } return $r; }
sub supt { my $a = shift; my $r = szero(); for my $n (1..$NT) { $r->[$n] = $a->[$n-1]->copy; } return $r; }      # t a
sub ssqrt1 { my $a = shift; die "ssqrt1: constant term must be 1" unless $a->[0] == 1; my $r = szero(); $r->[0] = rat(1);
  for my $n (1..$NT) { my $s = $a->[$n]->copy; for my $k (1..$n-1) { $s -= $r->[$k] * $r->[$n-$k]; } $r->[$n] = $s / 2; } return $r; }
sub slog1 { my $a = shift; die "slog1: constant term must be 1" unless $a->[0] == 1; return sint(smul(sder($a), sinv($a))); }
sub fact { my $n = shift; my $f = rat(1); $f *= $_ for 1..$n; return $f; }
sub dfact { my $n = shift; my $r = rat(1); for (my $i = $n; $i > 0; $i -= 2) { $r *= $i; } return $r; }        # n!!, (-1)!! = 1
sub spoly { my ($poly, $p) = @_; my $s = szero(); my $pp = cst(1); for my $co (@$poly) { $s = sadd($s, ssc($co, $pp)); $pp = smul($pp, $p); } return $s; }

my (%nC, %nB);
# chk(category, label, ok [, number of coefficients compared, number of mismatching coefficients]); seq() returns the last three
sub chk { my ($cat, $lab, $ok, $n, $nb) = @_; $n //= 1; $nb //= ($ok ? 0 : 1); $nC{$cat} += $n;
  if ($nb) { $nB{$cat} += $nb; print "  MISMATCH [$cat] $lab ($nb)\n" if $nB{$cat} <= 8; } }
sub seq { my ($a, $b, $upto) = @_; my $nb = 0; for my $n (0..$upto) { $nb++ unless $a->[$n] == $b->[$n]; } return ($nb == 0, $upto + 1, $nb); }

# Gaussian-normalised formal Laplace expansion (S1 s.2.2, (D1)) to z^M, with series-valued data:
#   Lnorm[f, g; u0] = Sum_n (2n-1)!! (-z/f2)^n [v^(2n)] ( g(u0+v) Exp[R(v)/z] ),  R = Sum_{k>=3} f^(k) v^k/k!.
# $fk->{k} (k = 3..2M+2) are the derivatives f^(k)(u0) (higher ones cannot contribute to z^m, m <= M); $g->[a] = [v^a] g(u0+v).
sub lnorm { my ($M, $f2, $fk, $g) = @_;
  my $V = 6 * $M; my $J = 2 * $M; my $KM = 2 * $M + 2;
  my @Rv; for my $k (3..$KM) { $Rv[$k] = ssc(rat(1) / fact($k), $fk->{$k}); }
  my @E; $E[0] = []; $E[0][0] = cst(1);
  for my $j (1..$J) { my @e;
    for my $a (0..$V) { next unless defined $E[$j-1][$a];
      for my $k (3..$KM) { last if $a + $k > $V; my $p = smul($E[$j-1][$a], $Rv[$k]); $e[$a+$k] = defined $e[$a+$k] ? sadd($e[$a+$k], $p) : $p; } }
    for my $a (0..$V) { $e[$a] = ssc(rat(1) / $j, $e[$a]) if defined $e[$a]; }
    $E[$j] = \@e; }
  my $mi = sinv(ssc(-1, $f2));
  my @mip = (cst(1)); push @mip, smul($mip[-1], $mi) for 1..3*$M;
  my @out;
  for my $m (0..$M) { my $tot = szero();
    for my $j (0..2*$m) { my $n = $m + $j; my $deg = 2 * $n; my $cc = szero();
      for my $a (0..$deg) { next unless defined $g->[$a]; my $b = $deg - $a; next unless defined $E[$j][$b]; $cc = sadd($cc, smul($g->[$a], $E[$j][$b])); }
      $tot = sadd($tot, ssc(dfact(2*$n - 1), smul($mip[$n], $cc))); }
    push @out, $tot; }
  return \@out; }

# theta^n p as a polynomial in p, theta p = p^2 - p (p = (1 - aU)^-1):  P_0 = p, P_n = (p^2 - p) dP_{n-1}/dp
my @P = ([rat(0), rat(1)]);
for my $n (1..2*$MZ + 2) { my @d = map { $P[$n-1][$_] * $_ } 1..$#{$P[$n-1]}; my @r = map { rat(0) } 0..$#d + 2;
  for my $i (0..$#d) { $r[$i+2] += $d[$i]; $r[$i+1] -= $d[$i]; } push @P, \@r; }
# W_k = (U d/dU)^k W, W = x log(1 - QU) - log U + log(1 - U):  W_1 = x - x pQ - p1,  W_k = -x P_{k-1}(pQ) - P_{k-1}(p1) (k >= 2)
sub Wk_series { my ($x, $pQ, $p1, $kmax) = @_; my %W;
  $W{1} = ssub(ssub(cst($x), ssc($x, $pQ)), $p1);
  for my $k (2..$kmax) { $W{$k} = ssub(ssc(rat(0) - $x, spoly($P[$k-1], $pQ)), spoly($P[$k-1], $p1)); }
  return \%W; }

my @B2k = (undef, rat("1/6"), rat("-1/30"), rat("1/42"), rat("-1/30"), rat("5/66"));
sub hw_coef { my ($x, $T, $pow) = @_; my $c = $x + 1;          # z-coefficients 0..T of h_w(z)^pow, h_w = exp(Sum_k B_2k/(2k(2k-1)) z^(2k-1) (x^(1-2k) + 1 - c^(1-2k)))
  my @e = map { rat(0) } 0..$T;
  for my $k (1..5) { my $p = 2*$k - 1; last if $p > $T; $e[$p] = rat($pow) * $B2k[$k] / (2*$k * (2*$k - 1)) * ($x**(-$p) + 1 - $c**(-$p)); }
  my @r = (rat(1)); for my $n (1..$T) { my $s = rat(0); for my $k (1..$n) { $s += $e[$k] * $r[$n-$k] * $k; } push @r, $s / $n; }
  return \@r; }
sub pair { my ($x, $a0, $a1, $b0, $b1) = @_; my $c = $x + 1; return $x*$x * ($a0*$b1 + $a1*$b0) - $c*$x * $a1*$b1; }   # <a, b>

my $t0 = time;
# ======================= Part A: fixed t != 0 =======================
{ local $NT = 0;
  for my $pt (["2", "5/8", "2"], ["1/3", "7/12", "1"], ["-3", "-1/4", "-1"]) {
    my ($x, $tau, $delta) = map { rat($_) } @$pt; my $c = $x + 1; my $M = 5;
    die "delta" unless $delta * $delta == $c * $c - 4 * $x * $tau;
    my $Q = rat(1) / (1 - $tau); my $sst = 2 / ($c + $delta); my $s2 = ($c + $delta) / (2 * $x * $tau);
    my $U1 = 1 - $tau * $sst; my $U2 = 1 - $tau * $s2;
    my $lab = "x=$x tau=$tau";
    chk("A1 roots", "$lab U2", $Q * $x * $U2**2 + $Q * (1 - $x) * $U2 - 1 == 0 && $U2 == ($x - 1 - $delta) / (2 * $x));
    chk("A1 roots", "$lab U1", $Q * $x * $U1**2 + $Q * (1 - $x) * $U1 - 1 == 0 && $U1 == ($x - 1 + $delta) / (2 * $x));
    chk("A1 roots", "$lab s2 critical", $x * $tau * $s2**2 - $c * $s2 + 1 == 0 && $s2 != 0 && $s2 != 1 && $s2 * $tau != 1);
    my $pQ = cst(1 / (1 - $Q * $U2)); my $p1 = cst(1 / (1 - $U2));
    my $W = Wk_series($x, $pQ, $p1, 2*$M + 2);
    chk("A1 roots", "$lab W_1(U_2) = 0", $W->{1}[0]->is_zero);
    my %fy = map { $_ => $W->{$_} } 3..2*$M+2;
    my $Ly = lnorm($M, $W->{2}, \%fy, [cst(1)]);
    # s-side: phi_tau^(k)(s) = (k-1)! [(-1)^(k-1) s^-k - x (1-s)^-k + tau^k (1 - tau s)^-k], amplitude (1 - tau s)^-1
    my %fs; my $f2s;
    for my $k (2..2*$M+2) { my $v = fact($k-1) * ((-1)**($k-1) * $s2**(-$k) - $x * (1 - $s2)**(-$k) + $tau**$k * (1 - $tau*$s2)**(-$k));
      if ($k == 2) { $f2s = $v; } else { $fs{$k} = cst($v); } }
    my @g = map { cst($tau**$_ * (1 - $tau*$s2)**(-$_ - 1)) } 0..6*$M;
    my $Ls = lnorm($M, cst($f2s), \%fs, \@g);
    for my $m (0..$M) { chk("A2 Step A at s_2 (z^0..z^5)", "$lab z^$m: " . $Ly->[$m][0] . " vs " . ($U2 * $Ls->[$m][0]), $Ly->[$m][0] == $U2 * $Ls->[$m][0]); }
    print "  A2 $lab: U_2 = $U2, s_2 = $s2, Wick row at U_2 = ", join(", ", map { $_->[0] } @$Ly), "\n";
    my $W2 = $W->{2}[0];
    my $L1 = ($x * $U1 + 1) / $tau; my $L2 = ($x * $U2 + 1) / $tau;
    my $D2q = $x**3 * ($L2 - $L1)**2 / ($L1 * ($c * $L1 - 2 * $x));      # 1/eta(e_2, e_2), e_2 = (H - L_1)/(L_2 - L_1)
    chk("A3 V1 at U_2, tensor rule", "$lab W2 = U2^2 phi''/tau^2", $W2 == $U2**2 * $f2s / $tau**2);
    chk("A3 V1 at U_2, tensor rule", "$lab W2 = x delta s*^2", $W2 == $x * $delta * $sst**2);
    chk("A3 V1 at U_2, tensor rule", "$lab phi''(s2) = delta/(x s2^2 U2^2)", $f2s == $delta / ($x * $s2**2 * $U2**2));
    chk("A3 V1 at U_2, tensor rule", "$lab Delta_2 = x^2 W2 = x^3 delta s*^2", $D2q == $x**2 * $W2 && $D2q == $x**3 * $delta * $sst**2);
    chk("A4 d f0/dt = L_2 - mu", "$lab", (1 - $tau) * $s2 / (1 - $tau * $s2) == 1 - $x * $s2 && 1 - $x * $s2 == $L2 - $x - $c * (1 - $tau) / $tau && $L2 == $x * $sst);
  }
}
printf "  part A done (%d s)\n", time - $t0;

# ======================= Part B: t-series =======================
for my $xs (@XS) {
  my $x = rat($xs); my $c = $x + 1; print "==================== x = $x (c = $c)\n";
  my $tau = szero(); for my $n (1..$NT) { $tau->[$n] = rat((-1)**($n+1)) / fact($n); }
  my $taut = szero(); for my $n (0..$NT) { $taut->[$n] = rat((-1)**$n) / fact($n+1); }               # tau/t
  my $ttau = sinv($taut);                                                                              # t/tau
  my $Qs = szero(); for my $n (0..$NT) { $Qs->[$n] = rat(1) / fact($n); }
  my $u = ssc(4 * $x / ($c * $c), $tau);
  my $d1 = ssqrt1(ssub(cst(1), $u)); my $delta = ssc($c, $d1);                                         # delta(0) = c
  my $sst = ssc(2, sinv(sadd(cst($c), $delta)));
  my $U1 = ssub(cst(1), smul($tau, $sst));
  my $U2 = ssc(rat(1) / (2 * $x), ssub(cst($x - 1), $delta));
  for my $pr (["U1", $U1], ["U2", $U2]) { my $U = $pr->[1];
    my $e = ssub(sadd(smul($Qs, ssc($x, smul($U, $U))), smul($Qs, ssc(1 - $x, $U))), cst(1));
    chk("B1 roots, L_2, Delta_2", "root $pr->[0]", seq($e, szero(), $NT)); }
  my $tL1 = smul(sadd(ssc($x, $U1), cst(1)), $ttau);                                                  # t L_1, L = (xU + 1)/tau
  my $tL2 = smul(sadd(ssc($x, $U2), cst(1)), $ttau);
  my $L2 = sdownt($tL2);                                                                               # exact to t^(NT-1)
  chk("B1 roots, L_2, Delta_2", "L_2 = x s*", seq($L2, ssc($x, $sst), $NT - 1));
  my $tt = supt(cst(1));                                                                               # the series t
  my $D2 = sdiv(ssc($x**3, smul(ssub($tL2, $tL1), ssub($tL2, $tL1))), smul($tL1, ssub(ssc($c, $tL1), ssc(2 * $x, $tt))));
  chk("B1 roots, L_2, Delta_2", "Delta_2 (quadratic form) = x^3 delta s*^2", seq($D2, ssc($x**3, smul($delta, smul($sst, $sst))), $NT));
  chk("B1 roots, L_2, Delta_2", "Delta_2(0) = x^3/c, Delta_2'(0) = 0", $D2->[0] == $x**3 / $c && $D2->[1]->is_zero);
  my $pQ = sinv(ssub(cst(1), smul($Qs, $U2))); my $p1 = sinv(ssub(cst(1), $U2));
  my $W = Wk_series($x, $pQ, $p1, 2*$MZ + 2);
  chk("B1 roots, L_2, Delta_2", "W_1(U_2) = 0", seq($W->{1}, szero(), $NT));
  chk("B1 roots, L_2, Delta_2", "V1 at U_2: Delta_2 = x^2 W_2(U_2)", seq($D2, ssc($x * $x, $W->{2}), $NT));
  # f2/tau^2 with f2 = phi''_tau(s_2), s_2 = 1/(x tau s*):  -x^2 s*^2 - x^3 s*^2/(1 - x tau s*)^2 + 1/U_2^2
  my $sst2 = smul($sst, $sst); my $om = ssub(cst(1), ssc($x, smul($tau, $sst)));
  my $f2t = sadd(ssub(ssc(rat(0) - $x*$x, $sst2), ssc($x**3, sdiv($sst2, smul($om, $om)))), sinv(smul($U2, $U2)));
  chk("B1 roots, L_2, Delta_2", "Delta_2 = x^2 U_2^2 phi''_tau(s_2)/tau^2", seq($D2, ssc($x*$x, smul(smul($U2, $U2), $f2t)), $NT));
  # B2: G := -c log(c s*) + x log(1 - x tau s*) - log(-x U_2) = phi_tau(s_2) + c log tau + const;  dG/dt = L_2 - x
  my $G = sadd(ssc(rat(0) - $c, slog1(ssc($c, $sst))), ssub(ssc($x, slog1($om)), slog1(ssc(rat(0) - $x, $U2))));
  chk("B2 d/dt[phi_tau(s_2) + c log tau] = L_2 - x", "", seq(sder($G), ssub($L2, cst($x)), $NT - 2));
  # rows at U_2 from the code's definitions
  my %fy = map { $_ => $W->{$_} } 3..2*$MZ+2;
  my $rpt = lnorm($MZ, $W->{2}, \%fy, [cst(1)]);                                                      # point frame, wickZ at U_2(t)
  my $hinv = hw_coef($x, $MZ, -1); my $hinv2 = hw_coef($x, $MZ, -2); my $H2 = hw_coef($x, $MZ, 2);
  my @r0; for my $m (0..$MZ) { my $s = szero(); for my $j (0..$m) { $s = sadd($s, ssc($hinv->[$m-$j], $rpt->[$j])); } $r0[$m] = $s; }
  my $dl = sdiv(sder($D2), ssc(2, $D2));                                                               # Delta_2'/(2 Delta_2)
  my @r1; $r1[0] = smul($L2, $r0[0]);
  for my $m (1..$MZ) { $r1[$m] = sadd(smul($L2, $r0[$m]), ssub(sder($r0[$m-1]), smul($r0[$m-1], $dl))); }
  for my $m (1..$MZ) {
    my $lhs = supt(ssub(sder($r1[$m-1]), smul($r1[$m-1], $dl)));
    my $rhs = smul($tL1, ssub($r1[$m], smul($L2, $r0[$m])));
    chk("B3 Row 2 at U_2 (Lemma 10.5), z^1..z^$MZ, t^0..t^$MAXC", "m = $m", seq($lhs, $rhs, $MAXC)); }
  for my $m (0..$MZ) {
    chk("B4 10.6(a): r^pt(z;0) = h^-1, r(z;0) = h^-2", "m = $m", $rpt->[$m][0] == $hinv->[$m] && $r0[$m][0] == $hinv2->[$m]);
    chk("B4 10.6(b) and [t^1] r_{0,2} = 0", "m = $m", $r1[$m][0] == $x / $c * $r0[$m][0] && $r0[$m][1]->is_zero); }
  # E_22 at t = 0 from (1.1): E(m, n) = Sum_{k=0}^{n} (-1)^k P_{m+1+k, n-k},  P_pq = -(-1)^(p+q) <v_2[p], v_2[q]>
  my @v0 = map { $r0[$_][0] } 0..$MZ; my @v1 = map { $r1[$_][0] } 0..$MZ;
  my $Ppq = sub { my ($p, $q) = @_; return rat(0) - rat((-1)**($p+$q)) * pair($x, $v0[$p], $v1[$p], $v0[$q], $v1[$q]); };
  my $Npq = sub { my ($p, $q) = @_; return (($p == 0 && $q == 0) ? rat(1) : rat(0)) - $H2->[$p] * $H2->[$q]; };   # [z^p w^q] (1 - H(z) H(w))
  my (%E22, $agreeNote, $nEnt); $agreeNote = 0; $nEnt = 0;
  for my $m (0..$MZ-1) { for my $n (0..$MZ-1-$m) {
    my $e = rat(0); my $f = rat(0);
    for my $k (0..$n) { $e += rat((-1)**$k) * $Ppq->($m+1+$k, $n-$k); $f += rat((-1)**$k) * $Npq->($m+1+$k, $n-$k); }
    $f *= $x**3 / $c; $E22{"$m,$n"} = $e; $nEnt++;
    chk("B5 E_22(t=0) = (x^3/c)(1 - h^2 h^2)/(z+w) (10.6(c))", "($m,$n)", $e == $f);
    $agreeNote++ if $e == (rat(0) - $x**3 / 2) / $c**($m+$n+2); } }
  chk("B5 E_22(0,0), minor", "E22(0,0) = -x^2(1+x+x^2)/(6c^2)", $E22{"0,0"} == (rat(0) - $x*$x * (1 + $x + $x*$x)) / (6 * $c*$c));
  my $S1 = $x**(-1) + 1 - $c**(-1); my $S2 = $x**(-3) + 1 - $c**(-3);
  my $minor = $E22{"0,0"} * $E22{"1,1"} - $E22{"0,1"} * $E22{"1,0"};
  chk("B5 E_22(0,0), minor", "minor = (x^6/c^2) S1 (72 S2 + 5 S1^3)/77760 != 0", $minor == $x**6 / $c**2 * $S1 * (72 * $S2 + 5 * $S1**3) / 77760 && !$minor->is_zero);
  print "  E_22(0,0)|t=0 = ", $E22{"0,0"}, ", 2x2 minor = $minor; the note's F3 value -(x^3/2)/c^(m+n+2) agrees at $agreeNote of $nEnt entries\n";
  # B7: K_2 = <v_2(z), v_2(-z)>/Delta_2 = 1 (Remark 10.6''): [z^n] Sum_{p+q=n} (-1)^q <v_2[p], v_2[q]> = Delta_2 delta_{n0}, as t-series
  for my $n (0..$MZ) { my $s = szero();
    for my $p (0..$n) { my $q = $n - $p; my $pr = sadd(ssc($x*$x, sadd(smul($r0[$p], $r1[$q]), smul($r1[$p], $r0[$q]))), ssc(rat(0) - $c*$x, smul($r1[$p], $r1[$q])));
      $s = sadd($s, ssc(rat((-1)**$q), $pr)); }
    chk("B7 K_2 = 1 (Remark 10.6''), z^0..z^$MZ, t^0..t^$MAXC", "z^$n", seq($s, ($n == 0 ? $D2 : szero()), $MAXC)); }
  # B6: rho from the quadratic form of Delta_1 versus the branch-free closed form
  my $t2D1 = sdiv(ssc($x**3, smul(ssub($tL1, $tL2), ssub($tL1, $tL2))), smul($L2, ssub(ssc($c, $L2), cst(2 * $x))));   # exact to t^(NT-1)
  my $y = ssc(rat(1) / (rat(0) - $x * $c**3), $t2D1);
  chk("B6 CF.2(d) rho, branch-free form", "t^2 Delta_1/(-x c^3) = 1 at t = 0", $y->[0] == 1);
  local $NT = $NT - 1;
  my @yy = @{$y}[0..$NT]; my $rho = sinv(ssqrt1(\@yy));
  my @uu = @{$u}[0..$NT]; my @tt2 = @{$taut}[0..$NT];
  my $q12 = ssqrt1(ssub(cst(1), \@uu)); my $q14 = ssqrt1($q12);
  my $rhocf = ssc(2, smul(\@tt2, sinv(smul($q14, sadd(cst(1), $q12)))));
  chk("B6 CF.2(d) rho, branch-free form", "rho to t^$MAXC", seq($rho, $rhocf, $MAXC));
  print "  rho = ", join(", ", map { $rho->[$_] } 0..3), ", ...\n";
  printf "  done x = %s (%d s)\n", $x, time - $t0;
}
print "==================== summary (comparisons / mismatches)\n";
my $tot = 0; my $bad = 0;
for my $k (sort keys %nC) { printf "  %-66s %5d %3d\n", $k, $nC{$k}, $nB{$k} // 0; $tot += $nC{$k}; $bad += $nB{$k} // 0; }
print "  total comparisons: $tot (part A at three rational points with t != 0; part B at x = @XS)\n";
print $bad ? "  MISMATCHES: $bad\n" : "  all comparisons exact\n";
