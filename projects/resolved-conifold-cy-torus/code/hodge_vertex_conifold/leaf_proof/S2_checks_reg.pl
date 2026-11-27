#!/usr/bin/perl
# S2_checks_reg.pl -- exact rational spot checks (Perl + Math::BigRat, no Mathematica) for section S2.10 of the leaf proof
# (leaf_proof/S2_metric_qde.md): the regular idempotent U_2. No closed form of S1 is used anywhere in this script.
# Part P: exact jets in e = tau - tau0 at points tau0 != 0 where everything is rational (parametrised by the rational regular
#   root s~(tau0) of x tau s^2 - c s + 1):
#   - the code's rows at both critical points, built as rc_lib.m / hvc_lib.m build them: Wick series (wickDirect) in y = log U with
#     W_k = (U d/dU)^k W, frame-B twist by h_w(z)^-1, literal Delta_i = 1/(e_i^T eta e_i), H-row with d/dt = (1 - tau) d/dtau;
#   - Lr := Lnorm[phi_tau, (1 - tau s)^-1; s~], the Wick formula (D1) of S1 in the coordinate s with amplitude (1 - tau s)^-1;
#   checks: the data of Lemma S2.10.1; (A2) r^pt_{0,2} = U_2 Lr (Lemma S2.10.2); (B2) the Gauss identity
#   [tau (1 - tau) Dop^2 + (gamma - (alpha + beta + 1) tau) Dop - alpha beta] Lr = 0, Dop = d/dtau + f0_tau/z - f2_tau/(2 f2)
#   (Lemma S2.10.3), coefficientwise in z from z^-2; the second QDE row at both idempotents ((V5); Props S2.7, S2.10.4);
#   <v_i(w), v_j(-w)> = delta_ij Delta_i ((V6); Props S2.9, S2.10, S2.10.5). Negative control: the sign of f2_tau/(2 f2) flipped.
# Part T: t-series at t = 0 for the regular idempotent (every ingredient is analytic there): L_2 = x s*, Delta_2 (literal) =
#   x^2 W2(U_2) = x^3 delta s*^2 with Delta_2(0) = x^3/c and Delta_2'(0) = 0; the second QDE row at U_2; v_2(z; 0) = h_w(z)^-2 (1, x/c)
#   (S4 Prop. 10.6(a), (b)); K_2 = 1, i.e. <v_2(w), v_2(-w)> = Delta_2; f2 = phi''_tau(s~) = Delta_2 tau^2/(x^2 U_2^2) =
#   (x^3/c) tau^2 (1 + O(tau)). Negative control: the sign of Delta_2'/(2 Delta_2) in the H-row flipped.
# Usage: perl S2_checks_reg.pl [M] [J] [N]   (M = max z-order, J = jet order, N = t-order; defaults 6 4 8)
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $M = @ARGV ? $ARGV[0] : 6;
my $J = @ARGV > 1 ? $ARGV[1] : 4;
my $N = @ARGV > 2 ? $ARGV[2] : 8;
my ($ntest, $nfail) = (0, 0); my %count;
sub ok { my ($cond, $msg) = @_; $ntest++; $count{(split /:/, $msg)[0]}++; unless ($cond) { $nfail++; print "FAIL: $msg\n"; } }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
sub rpow { my ($a, $n) = @_; my $r = R(1); if ($n >= 0) { $r *= $a for 1..$n; } else { $r /= $a for 1..(-$n); } $r }

# ---------------------------------------------------------------- truncated power series in one variable (indices 0..$LEN)
our $LEN;
sub sz { [ map { R(0) } 0..$LEN ] }
sub sc { my $a = sz(); $a->[0] = R($_[0]); $a }
sub sadd { my ($a, $b) = @_; [ map { $a->[$_] + $b->[$_] } 0..$LEN ] }
sub ssub { my ($a, $b) = @_; [ map { $a->[$_] - $b->[$_] } 0..$LEN ] }
sub ssc  { my ($s, $a) = @_; [ map { $a->[$_] * $s } 0..$LEN ] }
sub smul { my ($a, $b) = @_; my $c = sz();
  for my $i (0..$LEN) { next if $a->[$i]->is_zero; for my $k (0..$LEN-$i) { $c->[$i+$k] += $a->[$i] * $b->[$k]; } } $c }
sub sinv { my $a = shift; die "sinv: zero constant term" if $a->[0]->is_zero; my $b = sz(); $b->[0] = 1 / $a->[0];
  for my $n (1..$LEN) { my $s = R(0); $s += $a->[$_] * $b->[$n-$_] for 1..$n; $b->[$n] = -$s / $a->[0]; } $b }
sub sder { my $a = shift; my $b = sz(); $b->[$_-1] = $a->[$_] * $_ for 1..$LEN; $b }          # top coefficient becomes invalid
sub sdown { my $a = shift; die "sdown: nonzero constant term" unless $a->[0]->is_zero;           # divide by the variable
  my $b = sz(); $b->[$_-1] = $a->[$_] for 1..$LEN; $b }                                          # (top coefficient invalid)
sub sup { my $a = shift; my $b = sz(); $b->[$_] = $a->[$_-1] for 1..$LEN; $b }                   # multiply by the variable
sub ssqrt1 { my $a = shift; die "ssqrt1: constant term != 1" unless $a->[0] == 1; my $b = sz(); $b->[0] = R(1);
  for my $n (1..$LEN) { my $s = R(0); $s += $b->[$_] * $b->[$n-$_] for 1..$n-1; $b->[$n] = ($a->[$n] - $s) / 2; } $b }
sub seq { my ($a, $b, $upto) = @_; for my $i (0..$upto) { return 0 unless $a->[$i] == $b->[$i]; } 1 }
sub spowl { my ($a, $n) = @_; my @p = (sc(1)); push @p, smul($p[-1], $a) for 1..$n; \@p }       # [a^0, ..., a^n]

# ---------------------------------------------------------------- Wick sums
my %PCACHE;
sub parts { my ($n, $k, $hi) = @_;      # partitions of n into exactly k parts in [3, hi], non-increasing
  return ([]) if $k == 0 && $n == 0; return () if $k == 0 || $n < 3*$k;
  my $key = "$n,$k,$hi"; return @{$PCACHE{$key}} if exists $PCACHE{$key};
  my @out; my $top = $n - 3*($k-1); $top = $hi if $hi < $top;
  for my $p (reverse 3..$top) { next if $p*$k < $n; for my $rest (parts($n - $p, $k - 1, $p)) { push @out, [$p, @$rest]; } }
  $PCACHE{$key} = \@out; @out }
# Lnorm[f, g; u0] (S1, (D1)): [z^m] = Sum_a ga[a] Sum_n (2n-1)!! (-1/f2)^n Sum_{partitions of 2n - a into n - m parts >= 3} prod r_k^m_k/m_k!
# with r_k = f^(k)(u0)/k!, ga[a] = g^(a)(u0)/a!. For g = 1 (ga = [1]) this is wickDirect of rc_lib.m with Wk[k] -> f^(k)(u0).
sub wick { my ($f2, $r, $ga, $MM) = @_;
  my $if2 = sinv($f2); my @ip = (sc(1)); push @ip, smul($ip[-1], $if2) for 1..3*$MM;
  my @L;
  for my $m (0..$MM) { my $tot = sz();
    for my $a (0..2*$m) { next unless defined $ga->[$a];
      for my $n ($m..(3*$m - $a)) { my $j = $n - $m; my $b = 2*$n - $a; next if $b < 0;
        my $dfact = fact(2*$n) / (fact($n) * rpow(R(2), $n)); $dfact = -$dfact if $n % 2;
        my $inner;
        if ($j == 0) { next unless $b == 0; $inner = sc(1); }
        else { $inner = sz(); for my $p (parts($b, $j, $b)) { my %t; $t{$_}++ for @$p; my $mon = sc(1); my $den = R(1);
            for my $k (keys %t) { for (1..$t{$k}) { $mon = smul($mon, $r->[$k]); } $den *= fact($t{$k}); }
            $inner = sadd($inner, ssc(1/$den, $mon)); } }
        $tot = sadd($tot, ssc($dfact, smul($ga->[$a], smul($ip[$n], $inner)))); } }
    push @L, $tot; }
  \@L }

# ---------------------------------------------------------------- h_w(z)^sgn, h_w = exp(Sum_k B_2k/(2k(2k-1)) z^(2k-1) (x^(1-2k) + 1 - c^(1-2k)))
my @bern = (undef, R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"));
sub hw_coeffs { my ($x, $sgn, $MM) = @_; my $c = 1 + $x; my @ex = map { R(0) } 0..$MM;
  for my $k (1..8) { my $d = 2*$k - 1; last if $d > $MM;
    $ex[$d] = $sgn * $bern[$k] / (2*$k*(2*$k-1)) * (rpow($x, 1-2*$k) + 1 - rpow($c, 1-2*$k)); }
  my @h = (R(1), map { R(0) } 1..$MM);
  for my $n (1..$MM) { my $s = R(0); $s += $ex[$_] * $h[$n-$_] * $_ for 1..$n; $h[$n] = $s / $n; }
  \@h }

# ---------------------------------------------------------------- theta-polynomials: W_k = -P_k(A) - x P_k(B), A = (1-U)^-1, B = (1-QU)^-1
my $KM = 2*$M + 2; my @P; $P[1] = [R(0), R(1)];
for my $k (2..$KM) { my $p = $P[$k-1]; my @q = map { R(0) } 0..$k;
  for my $n (1..$#$p) { next if $p->[$n]->is_zero; $q[$n+1] += $p->[$n] * $n; $q[$n] -= $p->[$n] * $n; } $P[$k] = \@q; }
# the code's identity row at a critical point U: point frame (wickDirect) and frame B (twist by h_w^-1); also W_2(U)
sub code_rows { my ($x, $U, $Q, $hinv) = @_;
  my $A = sinv(ssub(sc(1), $U)); my $B = sinv(ssub(sc(1), smul($Q, $U)));
  my $Ap = spowl($A, $KM); my $Bp = spowl($B, $KM);
  my @W; for my $k (2..$KM) { my $w = sz(); for my $n (1..$k) { next if $P[$k][$n]->is_zero;
      $w = ssub($w, ssc($P[$k][$n], sadd($Ap->[$n], ssc($x, $Bp->[$n])))); } $W[$k] = $w; }
  my @r; $r[$_] = ssc(1/fact($_), $W[$_]) for 3..$KM;
  my $rpt = wick($W[2], \@r, [sc(1)], $M);
  my @r0; for my $m (0..$M) { my $v = sz(); $v = sadd($v, ssc($hinv->[$m-$_], $rpt->[$_])) for 0..$m; push @r0, $v; }
  return ($rpt, \@r0, $W[2]);
}

# ================================================================= Part P: jets in e = tau - tau0
$LEN = $J;
for my $pt ([R(2), R("5/8"), R(2)], [R("3/5"), R("16/27"), R("15/4")], [R(-3), R("5/12"), R(2)], [R("1/3"), R("7/12"), R(6)]) {
  my ($x, $tau0, $sig) = @$pt; my $c = 1 + $x;
  ok($x*$tau0*$sig*$sig - $c*$sig + 1 == 0 && (2*$x*$tau0*$sig - $c)/$c > 0, "P: base point: s~(tau0) is the root of x tau s^2 - c s + 1 with 2 x tau s - c = delta, delta/c > 0 (x=$x)");
  my $tau = sc($tau0); $tau->[1] = R(1);
  my $omt = ssub(sc(1), $tau); my $Q = sinv($omt);                  # Q = 1/(1 - tau)
  my $s = sc($sig);                                                   # s~(tau) by Newton's method
  for (1..6) { my $F = sadd(ssub(ssc($x, smul($tau, smul($s, $s))), ssc($c, $s)), sc(1));
    my $Fp = ssub(ssc(2*$x, smul($tau, $s)), sc($c)); $s = ssub($s, smul($F, sinv($Fp))); }
  ok(seq(sadd(ssub(ssc($x, smul($tau, smul($s, $s))), ssc($c, $s)), sc(1)), sz(), $J), "P: x tau s~^2 - c s~ + 1 = 0 as a jet (x=$x)");
  my $sst = sinv(smul(ssc($x, $tau), $s));                            # s* = 1/(x tau s~)
  my $U2 = ssub(sc(1), smul($tau, $s)); my $U1 = ssub(sc(1), smul($tau, $sst));
  ok(seq(sadd($U1, $U2), sc(($x - 1)/$x), $J) && seq(smul($U1, $U2), ssc(-1/$x, sinv($Q)), $J), "P: Vieta U_1 + U_2 = (x-1)/x, U_1 U_2 = -1/(Q x) (x=$x)");
  my $hinv = hw_coeffs($x, -1, $M);
  my ($rpt2, $r02, $W22) = code_rows($x, $U2, $Q, $hinv);
  my ($rpt1, $r01, $W21) = code_rows($x, $U1, $Q, $hinv);
  my $LofU = sub { my $UU = shift; smul(ssc(-1, smul($Q, sadd(ssc($x, $UU), sc(1)))), sinv(ssub(sc(1), $Q))) };   # rc_lib LofU
  my %L = (1 => $LofU->($U1), 2 => $LofU->($U2));
  my %Del; my %W2 = (1 => $W21, 2 => $W22);
  for my $i (1, 2) { my ($Li, $Lj) = ($L{$i}, $L{3-$i});                  # code: Delta_i = 1/(e_i^T eta e_i)
    my $ec = smul(ssub(ssc($c/rpow($x, 3), smul($Lj, $Lj)), ssc(2/rpow($x, 2), $Lj)), sinv(smul(ssub($Li, $Lj), ssub($Li, $Lj))));
    $Del{$i} = sinv($ec);
    ok(seq($Del{$i}, ssc($x*$x, $W2{$i}), $J), "P: (V1) Delta_$i = x^2 W2(U_$i) (x=$x)"); }
  # ---- data of Lemma S2.10.1 at the regular critical point
  my $del = ssub(ssc(2*$x, smul($tau, $s)), sc($c));                   # delta = 2 x tau s~ - c at s~
  ok(seq(smul($del, $del), ssub(sc($c*$c), ssc(4*$x, $tau)), $J), "P: delta := 2 x tau s~ - c, delta^2 = c^2 - 4 x tau (x=$x)");
  ok(seq(ssc($x, $s), $L{1}, $J) && seq(ssc($x, $sst), $L{2}, $J), "P: x s~ = L_1 and x s* = L_2 (code LofU) (x=$x)");
  my $oms = ssub(sc(1), $s); my $omts = $U2;                            # 1 - s~, 1 - tau s~ = U_2
  my $isp = spowl(sinv($s), $KM); my $iomsp = spowl(sinv($oms), $KM); my $taup = spowl($tau, $KM); my $iomtsp = spowl(sinv($omts), $KM + 1);
  my @fk; for my $k (2..$KM) { my $sg = ($k % 2) ? 1 : -1;               # phi_tau^(k)(s~) = (k-1)! [(-1)^(k-1) s^-k - x (1-s)^-k + tau^k (1 - tau s)^-k]
    $fk[$k] = ssc(fact($k-1), sadd(ssub(ssc($sg, $isp->[$k]), ssc($x, $iomsp->[$k])), smul($taup->[$k], $iomtsp->[$k]))); }
  my $f2 = $fk[2];
  ok(seq(smul(smul($omts, $omts), $f2), smul($del, sinv(ssc($x, smul($s, $s)))), $J), "P: (1 - tau s~)^2 f2 = delta/(x s~^2) (x=$x)");
  ok(seq($f2, smul($Del{2}, smul(smul($tau, $tau), sinv(ssc($x*$x, smul($U2, $U2))))), $J), "P: f2 = Delta_2 tau^2/(x^2 U_2^2) (x=$x)");
  my $f0t = smul($s, sinv($omts));                                      # d/dtau [phi_tau(s~(tau))] = s~/(1 - tau s~)
  my $lam2 = ssub(sc(1), ssc($x, $s));
  ok(seq(smul($omt, $f0t), $lam2, $J), "P: lambda_2 = (1 - tau) s~/(1 - tau s~) = 1 - x s~ (x=$x)");
  my $mu = ssub(ssc($c, sinv($tau)), sc(1));
  ok(seq(sadd($mu, $lam2), $L{2}, $J), "P: mu + lambda_2 = L_2, mu = c/tau - 1 (x=$x)");
  ok(seq($mu, sadd(sc($x), ssc($c, smul($omt, sinv($tau)))), $J), "P: mu = x + c (1 - tau)/tau (x=$x)");
  # ---- Lr and (A2)
  my @rk; $rk[$_] = ssc(1/fact($_), $fk[$_]) for 3..$KM;
  my @ga; $ga[$_] = smul($taup->[$_], $iomtsp->[$_+1]) for 0..2*$M;      # g^(a)(s~)/a! = tau^a (1 - tau s~)^(-a-1)
  my $Lr = wick($f2, \@rk, \@ga, $M);
  for my $m (0..$M) { ok(seq($rpt2->[$m], smul($U2, $Lr->[$m]), $J), "P: (A2) r^pt_{0,2}[$m] = U_2 Lr[$m] (x=$x)"); }
  # ---- (B2): Gauss identity for Lr, coefficientwise in z for z^-2 .. z^(M-2)
  my $hf = smul(sder($f2), sinv(ssc(2, $f2)));                           # f2_tau/(2 f2)
  my %X; $X{$_} = $Lr->[$_] for 0..$M;
  my $Dop = sub { my ($Xh, $sgnh) = @_; my %Y;
    for my $m (-2..$M) { my $xm = $Xh->{$m} // sz(); my $xm1 = $Xh->{$m+1} // sz();
      $Y{$m} = sadd(ssub(sder($xm), ssc($sgnh, smul($hf, $xm))), smul($f0t, $xm1)); } \%Y };
  for my $sgnh (1, -1) {
    my $D1 = $Dop->(\%X, $sgnh); my $D2 = $Dop->($D1, $sgnh); my $allok = 1;
    for my $m (-2..$M-2) {
      my $g = smul(smul($tau, $omt), $D2->{$m});
      $g = sadd($g, smul(ssub(sc($c), ssc(2, $tau)), $D1->{$m+1}));
      $g = sadd($g, smul(ssub(sc(2), ssc(3, $tau)), $D1->{$m}));
      my $xm = $X{$m} // sz(); my $xm1 = $X{$m+1} // sz(); my $xm2 = $X{$m+2} // sz();
      $g = ssub($g, sadd($xm2, sadd(ssc(2, $xm1), $xm)));
      my $okm = seq($g, sz(), $J - 2);
      if ($sgnh == 1) { ok($okm, "P: (B2) [tau(1-tau) Dop^2 + (gamma - (alpha+beta+1) tau) Dop - alpha beta] Lr = 0 at z^$m (x=$x)"); }
      else { $allok &&= $okm; } }
    ok(!$allok, "P: negative control: (B2) fails with the sign of f2_tau/(2 f2) flipped (x=$x)") if $sgnh == -1; }
  # ---- the code's H-row (d/dt = (1 - tau) d/dtau), second QDE row (V5) and symplecticity (V6) at both idempotents
  my $dt = sub { smul($omt, sder($_[0])) };
  my %r0 = (1 => $r01, 2 => $r02); my %r1; my %hl;
  for my $i (1, 2) { $hl{$i} = smul($dt->($Del{$i}), sinv(ssc(2, $Del{$i})));
    my @a = (smul($L{$i}, $r0{$i}[0]));
    for my $m (1..$M) { push @a, sadd(smul($L{$i}, $r0{$i}[$m]), ssub($dt->($r0{$i}[$m-1]), smul($r0{$i}[$m-1], $hl{$i}))); }
    $r1{$i} = \@a; }
  for my $i (1, 2) { my ($Li, $Lj) = ($L{$i}, $L{3-$i});
    for my $m (1..$M) { my $lhs = ssub($dt->($r1{$i}[$m-1]), smul($r1{$i}[$m-1], $hl{$i}));
      my $rhs = smul($Lj, ssub($r1{$i}[$m], smul($Li, $r0{$i}[$m])));
      ok(seq($lhs, $rhs, $J - 2), "P: (V5) second QDE row at U_$i, z^$m (x=$x)"); } }
  my $pair = sub { my ($a0, $a1, $b0, $b1) = @_; ssub(ssc($x*$x, sadd(smul($a0, $b1), smul($a1, $b0))), ssc($c*$x, smul($a1, $b1))) };
  for my $i (1, 2) { for my $j (1, 2) { for my $n (0..$M) { my $sum = sz();
        for my $q (0..$n) { my $p = $n - $q; my $tt = $pair->($r0{$i}[$p], $r1{$i}[$p], $r0{$j}[$q], $r1{$j}[$q]);
          $sum = ($q % 2) ? ssub($sum, $tt) : sadd($sum, $tt); }
        my $target = ($n == 0 && $i == $j) ? $Del{$i} : sz();
        ok(seq($sum, $target, $J - 1), "P: (V6) <v_$i(w), v_$j(-w)> = delta_ij Delta_i at w^$n (x=$x)"); } } }
  print "Part P, x = $x, tau0 = $tau0, s~ = $sig: $ntest checks so far, $nfail failures\n";
}

# ================================================================= Part T: t-series at t = 0, regular idempotent
$LEN = $N + 6;
for my $x (R(2), R("3/5"), R(-3), R("1/3")) { my $c = 1 + $x;
  my $tt = sz(); $tt->[$_] = R((($_ % 2) ? -1 : 1)) / fact($_ + 1) for 0..$LEN;      # tau/t = (1 - e^-t)/t
  my $tau = sup($tt); my $omt = ssub(sc(1), $tau);
  my $Q = sz(); $Q->[$_] = 1/fact($_) for 0..$LEN;                                    # e^t
  my $dc = ssqrt1(ssub(sc(1), ssc(4*$x/($c*$c), $tau)));                              # delta/c
  my $sst = ssc(2/$c, sinv(sadd(sc(1), $dc)));                                        # s* = 2/(c + delta)
  my $U2 = ssub(sc(1), ssc($c/(2*$x), sadd(sc(1), $dc)));                              # U_2 = 1 - (c + delta)/(2x)
  ok($U2->[0] == -1/$x, "T: U_2(0) = -1/x (x=$x)");
  ok(seq(sadd(ssub(ssc($x, smul($Q, smul($U2, $U2))), ssc($x - 1, smul($Q, $U2))), sc(-1)), sz(), $N), "T: LeqU at U_2 (x=$x)");
  my $hinv = hw_coeffs($x, -1, $M);
  my ($rpt2, $r02, $W22) = code_rows($x, $U2, $Q, $hinv);
  my $L2 = smul(ssc(-1, $Q), smul(sdown(sadd(ssc($x, $U2), sc(1))), sinv(sdown(ssub(sc(1), $Q)))));   # LofU at U_2
  ok(seq($L2, ssc($x, $sst), $N) && $L2->[0] == $x/$c, "T: L_2 = x s*, L_2(0) = x/c (x=$x)");
  my $tL2 = smul($tau, $L2); my $tL1 = ssub(sc($c), $tL2);                              # tau L_1 = c - tau L_2 (Ssum = c/tau)
  ok(seq($tL1, ssc($x, ssub(sc(1), $U2)), $N), "T: tau L_1 = x (1 - U_2) (LofU at U_1 = (x-1)/x - U_2) (x=$x)");
  my $ec = smul(ssub(ssc($c/rpow($x, 3), smul($tL1, $tL1)), ssc(2/rpow($x, 2), smul($tau, $tL1))), sinv(smul(ssub($tL2, $tL1), ssub($tL2, $tL1))));
  my $D2 = sinv($ec);                                                                  # code: Delta_2 = 1/(e_2^T eta e_2)
  ok(seq($D2, ssc($x*$x, $W22), $N), "T: (V1) Delta_2 = x^2 W2(U_2) (x=$x)");
  ok(seq($D2, ssc(rpow($x, 3)*$c, smul($dc, smul($sst, $sst))), $N), "T: Delta_2 = x^3 delta s*^2 (x=$x)");
  ok($D2->[0] == rpow($x, 3)/$c && $D2->[1] == 0, "T: Delta_2(0) = x^3/c, Delta_2'(0) = 0 (x=$x)");
  ok(rpow($x, 2)*2*($x/$c) - $c*$x*($x/$c)**2 == rpow($x, 3)/$c, "T: <(1, x/c), (1, x/c)> = x^3/c (x=$x)");
  # f2 = phi''_tau(s~) with 1/s~ = x tau s*, 1/(1 - s~) = (1/s~)/(1/s~ - 1), tau/(1 - tau s~) = tau/U_2
  my $is = smul(ssc($x, $tau), $sst); my $ioms = smul($is, sinv(ssub($is, sc(1)))); my $tU = smul($tau, sinv($U2));
  my $f2 = sadd(ssub(ssc(-1, smul($is, $is)), ssc($x, smul($ioms, $ioms))), smul($tU, $tU));
  ok(seq($f2, smul($D2, smul(smul($tau, $tau), sinv(ssc($x*$x, smul($U2, $U2))))), $N), "T: f2 = Delta_2 tau^2/(x^2 U_2^2) (x=$x)");
  ok($f2->[0] == 0 && $f2->[1] == 0 && $f2->[2] == rpow($x, 3)/$c, "T: f2 = (x^3/c) t^2 + O(t^3), i.e. (x^3/c) tau^2 (1 + O(tau)) (x=$x)");
  # H-row and the second QDE row at U_2, multiplied by tau: tau (D_2 r_{1,2})[m-1] = (tau L_1)(r_{1,2}[m] - L_2 r_{0,2}[m])
  for my $sgnh (1, -1) {
    my $hl = ssc($sgnh, smul(sder($D2), sinv(ssc(2, $D2))));
    my @r1 = (smul($L2, $r02->[0]));
    for my $m (1..$M) { push @r1, sadd(smul($L2, $r02->[$m]), ssub(sder($r02->[$m-1]), smul($r02->[$m-1], $hl))); }
    my $allok = 1;
    for my $m (1..$M) { my $lhs = smul($tau, ssub(sder($r1[$m-1]), smul($r1[$m-1], $hl)));
      my $rhs = smul($tL1, ssub($r1[$m], smul($L2, $r02->[$m])));
      if ($sgnh == 1) { ok(seq($lhs, $rhs, $N), "T: (V5) second QDE row at U_2 to t^$N, z^$m (x=$x)"); }
      else { $allok &&= seq($lhs, $rhs, $N); } }
    if ($sgnh == -1) { ok(!$allok, "T: negative control: row 2 fails with the sign of Delta_2'/(2 Delta_2) flipped (x=$x)"); next; }
    # t = 0 values and K_2 = 1
    my $hinv2 = hw_coeffs($x, -2, $M);
    for my $m (0..$M) {
      ok($r02->[$m][0] == $hinv2->[$m], "T: r_{0,2}[$m](0) = [z^$m] h_w^-2 (S4 10.6(a)) (x=$x)");
      ok($r1[$m][0] == ($x/$c) * $hinv2->[$m], "T: r_{1,2}[$m](0) = (x/c) [z^$m] h_w^-2 (S4 10.6(b)) (x=$x)"); }
    my $pair = sub { my ($a0, $a1, $b0, $b1) = @_; ssub(ssc($x*$x, sadd(smul($a0, $b1), smul($a1, $b0))), ssc($c*$x, smul($a1, $b1))) };
    for my $n (0..$M) { my $sum = sz();
      for my $q (0..$n) { my $p = $n - $q; my $t2 = $pair->($r02->[$p], $r1[$p], $r02->[$q], $r1[$q]);
        $sum = ($q % 2) ? ssub($sum, $t2) : sadd($sum, $t2); }
      ok(seq($sum, ($n == 0 ? $D2 : sz()), $N), "T: K_2 = 1: <v_2(w), v_2(-w)> = Delta_2 at w^$n to t^$N (x=$x)"); }
  }
  print "Part T, x = $x: $ntest checks so far, $nfail failures\n";
}

print "S2_checks_reg.pl, M = $M, J = $J, N = $N: $ntest checks, $nfail failures\n";
print "  $_: $count{$_}\n" for sort keys %count;
