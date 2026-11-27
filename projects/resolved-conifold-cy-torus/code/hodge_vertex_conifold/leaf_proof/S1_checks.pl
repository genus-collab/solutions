#!/usr/bin/perl
# S1_checks.pl -- exact spot checks (Math::BigRat) for leaf_proof/S1_thimble.md.  No Mathematica.
# Conventions: eps1 = x, eps2 = 1, c = 1 + x; W(U) = x Log[1 - Q U] - Log[U] + Log[1 - U]; tau = 1 - 1/Q;
# phi_tau(s) = Log[s] + x Log[1 - s] - Log[1 - tau s]; s* = 2/(c + delta), delta = Sqrt[c^2 - 4 x tau].
# Lnorm[f, g; u0] = Sum_n (2n-1)!! (-z/f2)^n [v^(2n)] ( g(u0+v) Exp[R(v)/z] ),  R = f(u0+v) - f(u0) - f2 v^2/2.
# Checks:
#  C1  coordinate invariance at fixed t != 0: the code's Wick row (coordinate y = log U, amplitude 1, derivatives
#      W_k = (U d/dU)^k W at U1) equals N = (1 - tau s*) Lnorm[phi_tau, (1 - tau s)^-1; s*]  (z^0..z^M);
#  C2  A(z) := Lnorm[phi_0, 1; 1/c] equals h_w(z) c/(c+z), h_w = Exp[Sum_k B_2k/(2k(2k-1)) z^(2k-1) (x^(1-2k) + 1 - c^(1-2k))];
#  C3  N(z, tau) = A(z) F(z, tau) Exp[-dphi(tau)/z] alg(tau) in Q((z))[[tau]] mod (tau^K, z^(M+1)), and the negative
#      z-powers of the right-hand side vanish (regularity), N computed by the Wick formula in s with tau-series coefficients;
#      also the frame-B row h_w^-1 N at tau^0 and tau^1 against c/(c+z) and (1-x)^2 z/(c (c+z)(c+2z));
#  C4  the amplitude form of Euler's integration by parts (Gauss operator = total derivative) at 6 random rational points
#      (draws at which the identity is undefined are discarded and redrawn until 6 valid points have been tested);
#  C5  eikonal equation, P = delta, R = -x (c + 3 delta)/(delta (c + delta)), Phi_0 = 1/alg (transport equation) at a point;
#  C6  homogeneity (Lemma 1.1, Lemma 2.6, Prop. 3.1, Remark 3.2 of S1_thimble.md), as tau-series at x = 2, 1/3, -5/2:
#      w_k := tau^k W_k(U1) from (1.11) is analytic at tau = 0 with w_k(0) = (-1)^k phi_0^(k)(1/c) (k = 2..12), w_2 = -delta/(x s*^2);
#      in ytilde = -Log(1 - tau s)/tau the derivatives of phi_tau o sigma_tau at ytilde1 are (-1)^k w_k; the code's Wick row
#      wickZ(W_k(U1)) = wickZ(w_k) equals N(z, tau) mod (tau^7, z^5), tau = 0 included, with value A(z) = h_w c/(c+z) at tau = 0;
#      and Lnorm[phi_tau o sigma_tau, 1; ytilde1] = N(z, tau) (the uniform version of Prop. 3.1).
# Usage: perl S1_checks.pl
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
our $K = 1;                                  # truncation of tau-series (number of coefficients); set per check
sub zero { [ map { R(0) } 1..$K ] }
sub cst  { my $a = shift; my $s = zero(); $s->[0] = ref $a ? $a->copy : R($a); $s }
sub sadd { my ($a, $b) = @_; [ map { $a->[$_] + $b->[$_] } 0..$K-1 ] }
sub ssub { my ($a, $b) = @_; [ map { $a->[$_] - $b->[$_] } 0..$K-1 ] }
sub ssc  { my ($c, $a) = @_; [ map { $c * $a->[$_] } 0..$K-1 ] }
sub siz  { my $a = shift; for (@$a) { return 0 unless $_->is_zero } 1 }
sub smul { my ($a, $b) = @_; my $r = zero();
  for my $i (0..$K-1) { next if $a->[$i]->is_zero; for my $j (0..$K-1-$i) { $r->[$i+$j] += $a->[$i] * $b->[$j] } } $r }
sub sinv { my $a = shift; my $r = zero(); $r->[0] = 1 / $a->[0];
  for my $n (1..$K-1) { my $s = R(0); $s += $a->[$_] * $r->[$n-$_] for 1..$n; $r->[$n] = -$s / $a->[0] } $r }
sub stau { my $s = zero(); $s->[1] = R(1) if $K > 1; $s }
sub spow { my ($a, $n) = @_; my $r = cst(1); $r = smul($r, $a) for 1..$n; $r }
sub ssqrt1 { my $a = shift; die "ssqrt1: constant term must be 1" unless $a->[0] == 1; my $r = zero(); $r->[0] = R(1);
  for my $n (1..$K-1) { my $s = $a->[$n]->copy; $s -= $r->[$_] * $r->[$n-$_] for 1..$n-1; $r->[$n] = $s / 2 } $r }
sub sint { my $a = shift; my $r = zero(); $r->[$_] = $a->[$_-1] / $_ for 1..$K-1; $r }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
sub dfact { my $n = shift; my $r = R(1); for (my $i = $n; $i > 0; $i -= 2) { $r *= $i } $r }   # n!!, (-1)!! = 1

# Gaussian-normalised formal Laplace expansion to z^M.  $f2: series; $fk->{k} (k = 3..6M): series f^(k)(u0);
# $g->[a]: Taylor coefficient of v^a of the amplitude (series), a = 0..6M (missing = 0).  Returns [z^0 .. z^M].
sub lnorm { my ($M, $f2, $fk, $g) = @_;
  my $V = 6 * $M; my $J = 2 * $M;
  my @Rv = map { zero() } 0..$V;
  for my $k (3..$V) { $Rv[$k] = ssc(1 / fact($k), $fk->{$k}) }
  my @E; $E[0] = [ map { zero() } 0..$V ]; $E[0][0] = cst(1);
  for my $j (1..$J) { my @e = map { zero() } 0..$V;
    for my $a (0..$V) { next if siz($E[$j-1][$a]); for my $k (3..$V-$a) { $e[$a+$k] = sadd($e[$a+$k], smul($E[$j-1][$a], $Rv[$k])) } }
    $E[$j] = [ map { ssc(R(1) / $j, $_) } @e ] }
  my $mi = sinv(ssc(-1, $f2)); my @mip = (cst(1)); push @mip, smul($mip[-1], $mi) for 1..3*$M;
  my @out;
  for my $m (0..$M) { my $tot = zero();
    for my $j (0..2*$m) { my $n = $m + $j; my $deg = 2 * $n; my $cc = zero();
      for my $a (0..$deg) { next unless defined $g->[$a]; my $b = $deg - $a; next if $b > $V; $cc = sadd($cc, smul($g->[$a], $E[$j][$b])) }
      $tot = sadd($tot, ssc(dfact(2*$n-1), smul($mip[$n], $cc))) }
    push @out, $tot }
  \@out }

# phi_tau data at s* (series in tau, or numbers when K = 1): f_k, amplitude (1 - tau s)^-1, and N = (1 - tau s*) Lnorm
sub phidata { my ($x, $tau, $sst, $M) = @_;                 # $tau, $sst: series
  my $V = 6 * $M;
  my $sg = sinv($sst); my $om = sinv(ssub(cst(1), $sst)); my $nu = sinv(ssub(cst(1), smul($tau, $sst)));
  my %fk; my ($sgp, $omp, $nup, $tp) = (cst(1), cst(1), cst(1), cst(1));
  my @fall;
  for my $k (1..$V) { $sgp = smul($sgp, $sg); $omp = smul($omp, $om); $nup = smul($nup, $nu); $tp = smul($tp, $tau);
    my $v = ssc(fact($k-1), sadd(ssub(ssc((-1)**($k-1), $sgp), ssc($x, $omp)), smul($tp, $nup)));
    $fall[$k] = $v; $fk{$k} = $v if $k >= 3 }
  my @g; my $tpa = cst(1); my $nua = $nu;
  for my $a (0..$V) { $g[$a] = smul($tpa, $nua); $tpa = smul($tpa, $tau); $nua = smul($nua, $nu) }
  my $L = lnorm($M, $fall[2], \%fk, \@g);
  my $oms = ssub(cst(1), smul($tau, $sst));
  ([ map { smul($oms, $_) } @$L ], $fall[2], $fall[1]) }

my $bad = 0; sub chk { my ($ok, $msg) = @_; print(($ok ? "ok    " : "FAIL  "), $msg, "\n"); $bad++ unless $ok }

# ---------------- C1: coordinate invariance at fixed t != 0 ----------------
$K = 1;
for my $pt ([R(2), R("5/8"), R(2)], [R("1/3"), R("7/12"), R(1)]) {
  my ($x, $tau, $delta) = @$pt; my $c = 1 + $x; my $M = 5;
  die "delta" unless $delta * $delta == $c * $c - 4 * $x * $tau;
  my $sst = 2 / ($c + $delta); my $Q = 1 / (1 - $tau); my $U1 = 1 - $tau * $sst;
  chk($Q * $x * $U1**2 + $Q * (1 - $x) * $U1 - 1 == 0 && $U1 == ($x - 1 + $delta) / (2 * $x), "C1 x=$x tau=$tau: U1 = $U1 is the singular root of Q x U^2 + Q (1-x) U - 1");
  # y-coordinate: W_k = -x P_{k-1}(pQ) - P_{k-1}(p1), P_0 = p, P_n = (p^2 - p) P_{n-1}'
  my @P = ([R(0), R(1)]);
  for my $n (1..6*$M) { my @d = map { $P[$n-1][$_] * $_ } 1..$#{$P[$n-1]}; my @r = map { R(0) } 0..$#d + 2;
    for my $i (0..$#d) { $r[$i+2] += $d[$i]; $r[$i+1] -= $d[$i] } push @P, \@r }
  my $ev = sub { my ($poly, $p) = @_; my $s = R(0); my $pp = R(1); for my $co (@$poly) { $s += $co * $pp; $pp *= $p } $s };
  my $pQ = 1 / (1 - $Q * $U1); my $p1 = 1 / (1 - $U1);
  my %Wk; $Wk{$_} = cst(-$x * $ev->($P[$_-1], $pQ) - $ev->($P[$_-1], $p1)) for 2..6*$M;
  my $W2 = $Wk{2}; delete $Wk{2};
  my $Ly = lnorm($M, $W2, \%Wk, [cst(1)]);
  my ($Ns, $f2s) = phidata($x, cst($tau), cst($sst), $M);
  chk(1, "C1 W2(U1) = $W2->[0], (1 - tau s*)^2 phi''(s*)/tau^2 = " . ((1 - $tau * $sst)**2 * $f2s->[0] / $tau**2) . ", -delta/(x tau^2 s*^2) = " . (-$delta / ($x * $tau**2 * $sst**2)));
  chk($W2->[0] == (1 - $tau * $sst)**2 * $f2s->[0] / $tau**2 && $W2->[0] == -$delta / ($x * $tau**2 * $sst**2), "C1 W2(U1) = U1^2 phi''_tau(s*)/tau^2 = -delta/(x tau^2 s*^2)");
  my $okc = 1; for my $m (0..$M) { $okc = 0 unless $Ly->[$m][0] == $Ns->[$m][0] }
  chk($okc, "C1 x=$x tau=$tau: Wick row in y = log U equals N = (1 - tau s*) Lnorm[phi_tau,(1-tau s)^-1; s*], z^0..z^$M: " . join(", ", map { $_->[0] } @$Ly));
}

# ---------------- z-series helpers (reuse the series code with K = number of z-coefficients) ----------------
my @B2k = (undef, R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"), R("43867/798"));
sub sexp0 { my $a = shift; die "sexp0: constant term must be 0" unless $a->[0] == 0;   # exp of a series with zero constant term
  my $r = zero(); $r->[0] = R(1); for my $n (1..$K-1) { my $s = R(0); $s += $_ * $a->[$_] * $r->[$n-$_] for 1..$n; $r->[$n] = $s / $n } $r }
sub hw_series { my ($x, $T) = @_; my $c = 1 + $x; local $K = $T + 1; my $e = zero();     # h_w(z), z-coefficients 0..T
  for my $k (1..9) { my $p = 2*$k - 1; last if $p > $T; $e->[$p] = $B2k[$k] / (2*$k * (2*$k - 1)) * ($x**(-$p) + 1 - $c**(-$p)) }
  sexp0($e) }
sub geo_series { my ($c, $T) = @_; [ map { (-1)**$_ / $c**$_ } 0..$T ] }               # c/(c+z) = Sum (-z/c)^n

# ---------------- C2: A(z) = h_w(z) c/(c+z) ----------------
for my $x (R(2), R("1/3"), R(5), R("-3/7")) { my $c = 1 + $x; my $M = 7;
  $K = 1;
  my %fk; my @fall; for my $k (1..6*$M) { $fall[$k] = cst(fact($k-1) * $c**$k * ((-1)**($k-1) - $x**(1-$k))); $fk{$k} = $fall[$k] if $k >= 3 }
  my $A = lnorm($M, $fall[2], \%fk, [cst(1)]);
  chk($fall[2][0] == -$c**3 / $x, "C2 x=$x: phi_0''(1/c) = -c^3/x");
  my $hw = hw_series($x, $M); my $ge = geo_series($c, $M);
  my @st = map { my $m = $_; my $s = R(0); $s += $hw->[$_] * $ge->[$m-$_] for 0..$m; $s } 0..$M;
  my $ok = 1; for my $m (0..$M) { $ok = 0 unless $A->[$m][0] == $st[$m] }
  chk($ok, "C2 x=$x: Lnorm[phi_0, 1; 1/c] = h_w(z) c/(c+z) for z^0..z^$M (z^1: " . $A->[1][0] . ", z^2: " . $A->[2][0] . ")");
}

# ---------------- C3: N(z,tau) = A F Exp[-dphi/z] alg as tau-series; regularity; frame-B row ----------------
for my $x (R(2), R("1/3")) { my $c = 1 + $x; my $M = 4; my $KK = 5; my $T = $M + 2 * $KK;
  $K = $KK;
  my $tau = stau();
  my $u = ssc(-4 * $x / $c**2, $tau); my $dc = ssqrt1(sadd(cst(1), $u));                 # delta/c
  my $sst = sinv(ssc($c, sadd(cst(1), $dc))); $sst = ssc(2, $sst);                    # s* = 2/(c (1 + delta/c))
  my $crit = sadd(ssub(ssc($x, smul($tau, smul($sst, $sst))), ssc($c, $sst)), cst(1));
  chk(siz($crit), "C3 x=$x: x tau s*^2 - c s* + 1 = 0 (mod tau^$K)");
  my ($N, $f2, $f1) = phidata($x, $tau, $sst, $M);
  chk(siz($f1), "C3 x=$x: phi_tau'(s*) = 0 (mod tau^$K)");
  # closed form of phi''(s*): -delta/(x s*^2 (1 - tau s*)^2)
  my $oms = ssub(cst(1), smul($tau, $sst));
  my $f2cf = ssc(-$c / $x, smul($dc, sinv(smul(smul($sst, $sst), smul($oms, $oms)))));
  chk(siz(ssub($f2, $f2cf)), "C3 x=$x: phi''_tau(s*) = -delta/(x s*^2 (1 - tau s*)^2)");
  my $alg = smul($oms, ssqrt1(ssc(1 / $f2->[0], $f2)));
  my $algcf = smul(ssqrt1($dc), sinv(ssc($c, $sst)));                                 # alg = (delta/c)^(1/2)/(c s*)
  chk(siz(ssub($alg, $algcf)), "C3 x=$x: alg = (1 - tau s*) (phi''_tau(s*)/phi''_0(1/c))^(1/2) = (delta/c)^(1/2)/(c s*)");
  my $Sp = smul($sst, sinv($oms));                                                   # d(dphi)/dtau = s*/(1 - tau s*)
  my $Sp2 = smul(ssub(cst(1), ssc($x, $sst)), sinv(ssub(cst(1), $tau)));               # = (1 - x s*)/(1 - tau)
  chk(siz(ssub($Sp, $Sp2)), "C3 x=$x: s*/(1 - tau s*) = (1 - x s*)/(1 - tau)");
  my $dphi = sint($Sp);
  # Laurent-in-z arrays: hash zpow -> tau-series.  F, Exp[-dphi/z], A (to z^T), alg
  my %F; for my $n (0..$K-1) { local $K = $T + $n + 1; my $r = cst(1);
      for my $j (1..$n) { my $l = zero(); $l->[0] = R(1); $l->[1] = R($j) if $K > 1; $r = smul($r, smul($l, $l)) }
      for my $j (2..$n+1) { my $l = zero(); $l->[0] = $c; $l->[1] = R($j) if $K > 1; $r = smul($r, sinv($l)) }
      $r = ssc(1 / fact($n), $r);
      for my $p (0..$T + $n) { my $zp = $p - $n; my $ts = [ map { R(0) } 1..$KK ]; $ts->[$n] = $r->[$p]; $F{$zp} = defined $F{$zp} ? [ map { $F{$zp}[$_] + $ts->[$_] } 0..$KK-1 ] : $ts } }
  my %E; { my $pw = cst(1); for my $j (0..$K-1) { $E{-$j} = ssc((-1)**$j / fact($j), $pw); $pw = smul($pw, $dphi) } }
  my $hw = hw_series($x, $T); my $ge = geo_series($c, $T);
  my %A; for my $m (0..$T) { my $s = R(0); $s += $hw->[$_] * $ge->[$m-$_] for 0..$m; $A{$m} = cst($s) }
  my $lmul = sub { my ($X, $Y, $top) = @_; my %Z;
    for my $p (keys %$X) { for my $q (keys %$Y) { my $r = $p + $q; next if $r > $top; $Z{$r} = defined $Z{$r} ? sadd($Z{$r}, smul($X->{$p}, $Y->{$q})) : smul($X->{$p}, $Y->{$q}) } } \%Z };
  my $FE = $lmul->(\%F, \%E, $T); my $FEA = $lmul->($FE, \%A, $M); my %RHS = map { $_ => smul($FEA->{$_}, $alg) } keys %$FEA;
  my $okneg = 1; for my $p (-($K-1)..-1) { $okneg = 0 if defined $RHS{$p} && !siz($RHS{$p}) }
  chk($okneg, "C3 x=$x: [tau^k] (A F Exp[-dphi/z] alg) has no negative powers of z (k < $K)");
  my $okeq = 1; for my $m (0..$M) { $okeq = 0 unless siz(ssub($RHS{$m}, $N->[$m])) }
  chk($okeq, "C3 x=$x: N(z,tau) = A(z) F(z,tau) Exp[-dphi/z] alg(tau) mod (tau^$K, z^" . ($M+1) . ")");
  # frame B: h_w^-1 N; tau^0 -> c/(c+z); tau^1 -> (1-x)^2 z/(c (c+z)(c+2z))
  my $hinv; { local $K = $M + 1; $hinv = sinv(hw_series($x, $M)) }
  my @rB = map { my $m = $_; my $s = zero(); $s = sadd($s, ssc($hinv->[$_], $N->[$m-$_])) for 0..$m; $s } 0..$M;
  my $ok0 = 1; my $ok1 = 1;
  my $t1; { local $K = $M + 1; my $l1 = zero(); $l1->[0] = $c; $l1->[1] = R(1); my $l2 = zero(); $l2->[0] = $c; $l2->[1] = R(2);
    my $zz = zero(); $zz->[1] = R(1); $t1 = ssc((1 - $x)**2 / $c, smul($zz, sinv(smul($l1, $l2)))) }
  for my $m (0..$M) { $ok0 = 0 unless $rB[$m][0] == (-1)**$m / $c**$m; $ok1 = 0 unless $rB[$m][1] == $t1->[$m] }
  chk($ok0 && $ok1, "C3 x=$x: frame-B row h_w^-1 N: [tau^0] = c/(c+z), [tau^1] = (1-x)^2 z/(c (c+z)(c+2z)) to z^$M");
}

# ---------------- C4: Gauss operator on the Euler integrand = total s-derivative (amplitude form) ----------------
$K = 1;
srand(20261031);
{ my ($valid, $draws) = (0, 0);
  while ($valid < 6) { $draws++; die "C4: too many invalid draws" if $draws > 1000;
    my ($x, $s, $tau, $z) = map { R(int(rand(19)) - 9) / R(int(rand(7)) + 2) } 1..4;
    # discard draws where (4.1) is undefined (z = 0, s = 0, s = 1, tau s = 1), where tau = 0 or 1, or where it is trivially 0 = 0 (z = -1)
    next if $z == 0 || $z == -1 || $s == 0 || $s == 1 || $tau * $s == 1 || $tau == 0 || $tau == 1; my $c = 1 + $x;
    my $a = 1 / $z + 1; my $ga = $c / $z + 2; my $om = 1 - $tau * $s;
    my $lhs = $tau * (1 - $tau) * $a * ($a + 1) * $s**2 / $om**3 + ($ga - (2 * $a + 1) * $tau) * $a * $s / $om**2 - $a * $a / $om;
    my $k  = $s * (1 - $s) / $om**2; my $kp = (1 - 2 * $s) / $om**2 + 2 * $tau * $s * (1 - $s) / $om**3;
    my $php = 1 / $s - $x / (1 - $s) + $tau / $om;
    my $rhs = -$a * ($kp + $php * $k / $z);
    $valid++;
    chk($lhs == $rhs, "C4 point $valid: Gauss_z[e^(phi/z) (1-tau s)^-1] e^(-phi/z) = -(1/z+1)(k' + phi' k/z), k = s(1-s)(1-tau s)^-2, at (x,s,tau,z) = ($x,$s,$tau,$z)");
  }
  print "      C4: $valid valid points tested ($draws draws)\n";
}

# ---------------- C5: eikonal and transport at a point ----------------
{ my ($x, $tau, $delta) = (R(2), R("5/8"), R(2)); my $c = 1 + $x; my $s = 2 / ($c + $delta);
  my $sp = $x * $s**2 / $delta;                                     # ds*/dtau = x s*^2/delta
  my $Sp = (1 - $x * $s) / (1 - $tau); my $Spp = -$x * $sp / (1 - $tau) + (1 - $x * $s) / (1 - $tau)**2;
  my $eik = $tau * (1 - $tau) * $Sp**2 + ($c - 2 * $tau) * $Sp - 1;
  my $P = 2 * $tau * (1 - $tau) * $Sp + $c - 2 * $tau; my $Rr = $tau * (1 - $tau) * $Spp + (2 - 3 * $tau) * $Sp - 2;
  chk($eik == 0 && $P == $delta && $Rr == -$x * ($c + 3 * $delta) / ($delta * ($c + $delta)) && $sp == $x * $s**2 / ($c - 2 * $x * $tau * $s),
      "C5 x=2 tau=5/8: eikonal = 0, P = delta = $P, R = $Rr = -x(c+3delta)/(delta(c+delta)), s*' = x s*^2/delta");
  # Phi_0 = 1/alg solves P Phi_0' + R Phi_0 = 0: (log alg)' = R/P with alg = delta^(1/2) (c + delta)/(2 c^(3/2)):
  # (log alg)' = delta'/(2 delta) + delta'/(c + delta), delta' = -2x/delta
  my $dp = -2 * $x / $delta; my $lal = $dp / (2 * $delta) + $dp / ($c + $delta);
  chk($lal == $Rr / $P, "C5 x=2 tau=5/8: (log alg)' = R/P, i.e. Phi_0 = 1/alg solves the first transport equation");
}

# ---------------- C6: homogeneity, rescaled derivatives w_k = tau^k W_k(U1), the uniform coordinate ytilde ----------------
for my $x (R(2), R("1/3"), R("-5/2")) { my $c = 1 + $x; my $M = 4; my $KMAX = 12;   # KMAX >= 2M + 2: enough for z^0..z^M
  $K = 7;
  my $tau = stau();
  my $dc = ssqrt1(sadd(cst(1), ssc(-4 * $x / $c**2, $tau)));                          # delta/c, value 1 at tau = 0
  my $sst = ssc(2, sinv(ssc($c, sadd(cst(1), $dc))));                                   # s* = 2/(c + delta)
  my $U1 = ssub(cst(1), smul($tau, $sst));
  chk(siz(ssub(smul(sinv(ssub(cst(1), $tau)), sadd(ssc($x, smul($U1, $U1)), ssc(1 - $x, $U1))), cst(1))), "C6 x=$x: U1 = 1 - tau s* solves Q x U^2 + Q (1-x) U - 1 = 0 (Q = 1/(1-tau))");
  # (1.11): w_k = -Sum_i a_{k-1,i} tau^(k-i) [x A^i + B^i],  A := tau p_Q = -(1 - tau)/(1 - s*),  B := tau p_1 = 1/s*
  my $A = ssc(-1, smul(ssub(cst(1), $tau), sinv(ssub(cst(1), $sst))));
  my $B = sinv($sst);
  my @P = ([R(0), R(1)]);                                                              # P_0 = p, P_n = (p^2 - p) P_{n-1}'
  for my $n (1..$KMAX) { my @d = map { $P[$n-1][$_] * $_ } 1..$#{$P[$n-1]}; my @r = map { R(0) } 0..$#d + 2;
    for my $i (0..$#d) { $r[$i+2] += $d[$i]; $r[$i+1] -= $d[$i] } push @P, \@r }
  my %w;
  for my $k (2..$KMAX) { my $s = zero();
    die "P_{k-1}(0) != 0" unless $P[$k-1][0]->is_zero;
    for my $i (1..$#{$P[$k-1]}) { my $co = $P[$k-1][$i]; next if $co->is_zero; die "deg" if $i > $k;
      $s = ssub($s, ssc($co, smul(spow($tau, $k - $i), sadd(ssc($x, spow($A, $i)), spow($B, $i))))) }
    $w{$k} = $s }
  my $ok = 1; for my $k (2..$KMAX) { $ok = 0 unless $w{$k}[0] == R(-1)**$k * fact($k-1) * $c**$k * (R(-1)**($k-1) - $x**(1-$k)) }
  chk($ok, "C6 x=$x: w_k = tau^k W_k(U1) is analytic at tau = 0 with w_k(0) = (-1)^k phi_0^(k)(1/c), k = 2..$KMAX");
  chk(siz(ssub($w{2}, ssc(-$c / $x, smul($dc, sinv(smul($sst, $sst)))))), "C6 x=$x: w_2 = -delta/(x s*^2), i.e. (1.5)");
  # the y-side at tau != 0 is literally the code's W_k(U1): check (1.11) against (1.10), W_k = -x P_{k-1}(p_Q) - P_{k-1}(p_1), at (x, tau) = (2, 5/8)
  if ($x == 2) { my ($t0, $s0) = (R("5/8"), R("2/5")); my $Q0 = 1 / (1 - $t0); my $u0 = 1 - $t0 * $s0;
    my $ev = sub { my ($poly, $p) = @_; my $s = R(0); my $pp = R(1); for my $co (@$poly) { $s += $co * $pp; $pp *= $p } $s };
    my $okp = 1;
    for my $k (2..$KMAX) { my $Wk = -$x * $ev->($P[$k-1], 1 / (1 - $Q0 * $u0)) - $ev->($P[$k-1], 1 / (1 - $u0));
      my $wk = R(0); for my $i (1..$#{$P[$k-1]}) { $wk -= $P[$k-1][$i] * $t0**($k - $i) * ($x * (-(1 - $t0) / (1 - $s0))**$i + $s0**(-$i)) }
      $okp = 0 unless $wk == $t0**$k * $Wk }
    chk($okp, "C6 x=2 tau=5/8: (1.11) equals tau^k W_k(U1) computed from (1.10), W_k = -x P_{k-1}(1/(1-QU)) - P_{k-1}(1/(1-U)), k = 2..$KMAX") }
  # uniform coordinate ytilde = -Log(1 - tau s)/tau, s = sigma_tau(ytilde) = (1 - e^(-tau ytilde))/tau:
  # phi_tau(sigma_tau(ytilde1 + v)) = phi_tau(s* + Ds(v)),  Ds(v) = U1 (1 - e^(-tau v))/tau = U1 Sum_{j>=1} (-tau)^(j-1) v^j/j!
  my @Ds = (zero()); for my $j (1..$KMAX) { $Ds[$j] = ssc(R(1) / fact($j), smul($U1, spow(ssc(-1, $tau), $j - 1))) }
  my $sg = sinv($sst); my $om = sinv(ssub(cst(1), $sst)); my $nu = sinv($U1);
  my @phk; for my $k (1..$KMAX) { $phk[$k] = ssc(fact($k-1), sadd(ssub(ssc(R(-1)**($k-1), spow($sg, $k)), ssc($x, spow($om, $k))), smul(spow($tau, $k), spow($nu, $k)))) }
  my @f = map { zero() } 0..$KMAX; my @Dp = map { zero() } 0..$KMAX; $Dp[0] = cst(1);
  for my $k (1..$KMAX) { my @n = map { zero() } 0..$KMAX;
    for my $a (0..$KMAX) { next if siz($Dp[$a]); for my $b (1..$KMAX-$a) { $n[$a+$b] = sadd($n[$a+$b], smul($Dp[$a], $Ds[$b])) } }
    @Dp = @n;
    for my $d (0..$KMAX) { $f[$d] = sadd($f[$d], ssc(R(1) / fact($k), smul($phk[$k], $Dp[$d]))) } }
  $ok = siz($f[1]); for my $k (2..$KMAX) { $ok = 0 unless siz(ssub(ssc(fact($k), $f[$k]), ssc(R(-1)**$k, $w{$k}))) }
  chk($ok, "C6 x=$x: d^k/dytilde^k (phi_tau o sigma_tau)(ytilde1) = (-tau)^k W_k(U1) = (-1)^k w_k for k = 1..$KMAX (k = 1: 0)");
  # the code's Wick row in y = log U, rewritten by homogeneity (lambda = tau) as wickZ(w_k), against N(z, tau), tau = 0 included
  my %fw; for my $k (3..6*$M) { $fw{$k} = $k <= $KMAX ? $w{$k} : zero() }             # parts > 2M + 2 do not reach z^M
  my $Lw = lnorm($M, $w{2}, \%fw, [cst(1)]);
  my ($N) = phidata($x, $tau, $sst, $M);
  $ok = 1; for my $m (0..$M) { $ok = 0 unless siz(ssub($Lw->[$m], $N->[$m])) }
  chk($ok, "C6 x=$x: code's Wick row wickZ(W_k(U1)) = wickZ(w_k) equals N(z, tau) mod (tau^$K, z^" . ($M+1) . "), tau = 0 included");
  my $hw = hw_series($x, $M); my $ge = geo_series($c, $M);
  $ok = 1; for my $m (0..$M) { my $s = R(0); $s += $hw->[$_] * $ge->[$m-$_] for 0..$m; $ok = 0 unless $Lw->[$m][0] == $s }
  chk($ok, "C6 x=$x: r^pt_{0,1}(z; 0) = wickZ(w_k(0)) = A(z) = h_w(z) c/(c+z), z^0..z^$M");
  # uniform version of Prop. 3.1: Lnorm[phi_tau o sigma_tau, 1; ytilde1] = N(z, tau)
  my %fc; for my $k (3..6*$M) { $fc{$k} = $k <= $KMAX ? ssc(fact($k), $f[$k]) : zero() }
  my $Lc = lnorm($M, ssc(2, $f[2]), \%fc, [cst(1)]);
  $ok = 1; for my $m (0..$M) { $ok = 0 unless siz(ssub($Lc->[$m], $N->[$m])) }
  chk($ok, "C6 x=$x: uniform coordinate: Lnorm[phi_tau o sigma_tau, 1; ytilde1] = N(z, tau) mod (tau^$K, z^" . ($M+1) . ")");
}
print "\n", ($bad ? "$bad FAILURES" : "all checks passed"), "\n";
