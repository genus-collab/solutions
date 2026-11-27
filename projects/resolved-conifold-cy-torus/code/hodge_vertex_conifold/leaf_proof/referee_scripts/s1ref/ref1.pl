#!/usr/bin/perl
# Referee check R1-R4 for leaf_proof/S1_thimble.md (independent of S1_checks.pl).
# R1: the code's Wick row in the coordinate y = log U (wickDirect partition formula, amplitude 1), evaluated at the
#     singular root U1 as a TAU-SERIES (using homogeneity: W_k -> tau^k W_k leaves every Wick monomial unchanged),
#     equals A(z) F(z,tau) Exp[-dphi/z] alg(tau) mod (tau^K, z^(M+1)), with A = Exp[ex] c/(c+z) and ex built literally
#     as in hvc_lib.m hCoeffs: ex = -Sum_k B_2k/(2k(2k-1)) z^(2k-1) Sum_a w_a^(1-2k), w = (c, -x, -1).
# R2: frame-B row Exp[-ex] r^pt (hCoeffs with s = -1, as twistRows applies it) converted to t = log Q, then
#     rho q_t - q_0 against the stored leaf decompositions delta_1..delta_3 of leaf_structure.md (rho from leaf_structure.md).
# R3: rho(leaf_structure.md) = tau/(t alg) and Delta_1(leaf_structure.md Step 0) = -x delta/(tau^2 s*^2) (= x^2 W_2(U1)).
# R4: top coefficient a_{k+1}^(k) of leaf_structure.md Step 3 from the closed form, k = 1..K-2.
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
our $K;
sub z0 { [ map { R(0) } 1..$K ] }
sub cst { my $a = shift; my $s = z0(); $s->[0] = ref $a ? $a->copy : R($a); $s }
sub add { my ($a, $b) = @_; [ map { $a->[$_] + $b->[$_] } 0..$K-1 ] }
sub sbt { my ($a, $b) = @_; [ map { $a->[$_] - $b->[$_] } 0..$K-1 ] }
sub scl { my ($c, $a) = @_; [ map { $c * $a->[$_] } 0..$K-1 ] }
sub mul { my ($a, $b) = @_; my $r = z0();
  for my $i (0..$K-1) { next if $a->[$i]->is_zero; for my $j (0..$K-1-$i) { $r->[$i+$j] += $a->[$i] * $b->[$j] } } $r }
sub inv { my $a = shift; die "inv: zero constant" if $a->[0]->is_zero; my $r = z0(); $r->[0] = 1 / $a->[0];
  for my $n (1..$K-1) { my $s = R(0); $s += $a->[$_] * $r->[$n-$_] for 1..$n; $r->[$n] = -$s / $a->[0] } $r }
sub sqrt1 { my $a = shift; die "sqrt1" unless $a->[0] == 1; my $r = z0(); $r->[0] = R(1);
  for my $n (1..$K-1) { my $s = $a->[$n]->copy; $s -= $r->[$_] * $r->[$n-$_] for 1..$n-1; $r->[$n] = $s / 2 } $r }
sub integ { my $a = shift; my $r = z0(); $r->[$_] = $a->[$_-1] / $_ for 1..$K-1; $r }
sub pw { my ($a, $n) = @_; my $r = cst(1); $r = mul($r, $a) for 1..$n; $r }
sub iszero { my $a = shift; for (@$a) { return 0 unless $_->is_zero } 1 }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
# compose: f(tau) with tau = g(t), g(0) = 0
sub compose { my ($f, $g) = @_; my $r = z0(); my $gp = cst(1);
  for my $k (0..$K-1) { $r = add($r, scl($f->[$k], $gp)); $gp = mul($gp, $g) } $r }
# partitions of N into exactly J parts in [lo, hi], non-increasing
sub parts { my ($N, $J, $hi, $lo) = @_; return ($N == 0 ? ([]) : ()) if $J == 0; my @out;
  for my $p (reverse($lo..($hi < $N ? $hi : $N))) { next if $p * $J < $N; next if $lo * ($J-1) > $N - $p;
    for my $rest (parts($N - $p, $J - 1, $p, $lo)) { push @out, [$p, @$rest] } } @out }

my @B2k = (undef, R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"), R("43867/798"), R("-174611/330"));
my $bad = 0; sub chk { my ($ok, $msg) = @_; print(($ok ? "ok    " : "FAIL  "), $msg, "\n"); $bad++ unless $ok }

for my $x (R(2), R("1/3"), R("3/5"), R("-5/2")) {
  my $c = 1 + $x; my $M = 5; $K = 7; my $ZT = $M + $K + 1;
  my $T = z0(); $T->[1] = R(1);
  my $dc = sqrt1(add(cst(1), scl(-4 * $x / $c**2, $T)));            # delta/c
  my $sst = scl(2 / $c, inv(add(cst(1), $dc)));                      # s* = 2/(c + delta)
  chk(iszero(add(sbt(scl($x, mul($T, mul($sst, $sst))), scl($c, $sst)), cst(1))), "x=$x: x tau s*^2 - c s* + 1 = 0");
  my $U1 = sbt(cst(1), mul($T, $sst));
  # U1 must be the root (x - 1 + delta)/(2x) of Q x U^2 + Q (1-x) U - 1 (Q = 1/(1-tau)): check (1-tau)^-1 (x U1^2 + (1-x) U1) - 1 = 0
  my $Qs = inv(sbt(cst(1), $T));
  chk(iszero(sbt(mul($Qs, add(scl($x, mul($U1, $U1)), scl(1 - $x, $U1))), cst(1))) && iszero(sbt($U1, scl(1/(2*$x), add(cst($x - 1), scl($c, $dc))))), "x=$x: U1 = 1 - tau s* = (x-1+delta)/(2x) solves LeqU");
  # rescaled derivatives w_k = tau^k W_k(U1):  W_k = -x P_{k-1}(pQ) - P_{k-1}(p1), pQ = 1/(1 - Q U1), p1 = 1/(1 - U1)
  # tau pQ = -(1 - tau)/(1 - s*), tau p1 = 1/s*
  my $a = scl(-1, mul(sbt(cst(1), $T), inv(sbt(cst(1), $sst))));
  my $b = inv($sst);
  my @P = ([R(0), R(1)]);
  for my $n (1..2*$M+2) { my @d = map { $P[$n-1][$_] * $_ } 1..$#{$P[$n-1]}; my @r = map { R(0) } 0..$#d + 2;
    for my $i (0..$#d) { $r[$i+2] += $d[$i]; $r[$i+1] -= $d[$i] } push @P, \@r }
  my %w;
  for my $k (2..2*$M+2) { my $s = z0();
    for my $i (0..$#{$P[$k-1]}) { my $co = $P[$k-1][$i]; next if $co->is_zero; die "deg" if $i > $k;
      my $tp = pw($T, $k - $i);
      $s = sbt($s, scl($x * $co, mul(pw($a, $i), $tp)));
      $s = sbt($s, scl($co, mul(pw($b, $i), $tp))) }
    $w{$k} = $s }
  chk($w{2}[0] == -$c**3 / $x, "x=$x: tau^2 W_2(U1) -> -c^3/x = phi_0''(1/c)");
  # closed form (1.5): W_2(U1) = -delta/(x tau^2 s*^2)
  chk(iszero(sbt($w{2}, scl(-$c / $x, mul($dc, inv(mul($sst, $sst)))))), "x=$x: tau^2 W_2(U1) = -delta/(x s*^2)  [(1.5)]");
  my $w2i = inv($w{2});
  # wickDirect, literally: z^m coefficient = Sum_{n=m}^{3m} (2n)! (-1)^n/(n! 2^n W2^n) Sum_{partitions of 2n into n-m parts in [3, 2M+2]} prod W_k^m_k/((k!)^m_k m_k!)
  my @rpt;
  for my $m (0..$M) { my $tot = z0();
    for my $n ($m..3*$m) { my $j = $n - $m;
      my @pl = ($j == 0) ? ($n == 0 ? ([]) : ()) : parts(2*$n, $j, 2*$M+2, 3);
      for my $p (@pl) { my %mult; $mult{$_}++ for @$p;
        my $co = fact(2*$n) * R(-1)**$n / (fact($n) * R(2)**$n);
        my $mon = pw($w2i, $n);
        for my $k (keys %mult) { $co /= fact($k)**$mult{$k} * fact($mult{$k}); $mon = mul($mon, pw($w{$k}, $mult{$k})) }
        $tot = add($tot, scl($co, $mon)) } }
    push @rpt, $tot }
  # ex exactly as hCoeffs: p[m] = Total[w^-m], ex = -Sum_k B2k/(2k(2k-1)) z^(2k-1) p[2k-1]; z-series to z^ZT
  my @wts = ($c, -$x, R(-1));
  my @ex = map { R(0) } 0..$ZT;
  for my $k (1..10) { my $p = 2*$k - 1; last if $p > $ZT; my $pm = R(0); $pm += $_**(-$p) for @wts; $ex[$p] = -$B2k[$k] / (2*$k*(2*$k-1)) * $pm }
  my $zexp = sub { my ($e) = @_; my @r = (R(1)); for my $n (1..$ZT) { my $s = R(0); $s += $_ * $e->[$_] * $r[$n-$_] for 1..$n; $r[$n] = $s / $n } \@r };
  my $hw = $zexp->(\@ex);                       # S1's h_w = Exp[ex]
  my $hwinv = $zexp->([ map { -$_ } @ex ]);       # what twistRows multiplies by: Exp[s ex], s = -1
  # S1 section 0 form: h_w = Exp[Sum_k B2k/(2k(2k-1)) z^(2k-1) (x^(1-2k) + 1 - c^(1-2k))]
  my @ex2 = map { R(0) } 0..$ZT; for my $k (1..10) { my $p = 2*$k-1; last if $p > $ZT; $ex2[$p] = $B2k[$k]/(2*$k*(2*$k-1)) * ($x**(-$p) + 1 - $c**(-$p)) }
  my $ok = 1; for my $i (0..$ZT) { $ok = 0 unless $ex[$i] == $ex2[$i] } chk($ok, "x=$x: hCoeffs exponent ex equals Sum B2k/(2k(2k-1)) z^(2k-1) (x^(1-2k) + 1 - c^(1-2k))  [s.0]");
  my @A; for my $m (0..$ZT) { my $s = R(0); $s += $hw->[$_] * (R(-1)/$c)**($m-$_) for 0..$m; $A[$m] = $s }   # h_w c/(c+z)
  # RHS: Laurent in z (hash zpow -> tau-series)
  my %F;
  for my $n (0..$K-1) { my @S = (R(1)); push @S, R(0) for 1..$ZT+$n;
    my $pmul = sub { my ($q) = @_; my @r = map { R(0) } 0..$#S; for my $i (0..$#S) { next if $S[$i]->is_zero; for my $j (0..$#$q) { last if $i+$j > $#S; $r[$i+$j] += $S[$i] * $q->[$j] } } @S = @r };
    for my $j (1..$n) { $pmul->([R(1), R($j)]); $pmul->([R(1), R($j)]) }
    for my $j (2..$n+1) { $pmul->([ map { (R(-$j)/$c)**$_ / $c } 0..$#S ]) }
    for my $i (0..$#S) { my $zp = $i - $n; next if $zp > $ZT; my $ts = z0(); $ts->[$n] = $S[$i] / fact($n);
      $F{$zp} = defined $F{$zp} ? add($F{$zp}, $ts) : $ts } }
  my $dphi = integ(mul($sst, inv(sbt(cst(1), mul($T, $sst)))));       # d(dphi)/dtau = s*/(1 - tau s*)
  my %E; { my $q = cst(1); for my $l (0..$K-1) { $E{-$l} = scl(R(1)/fact($l), $q); $q = mul($q, scl(-1, $dphi)) } }
  my $alg = mul(sqrt1($dc), scl(R(1)/2, add(cst(1), $dc)));         # (delta/c)^(1/2) (c + delta)/(2c)
  my $lmul = sub { my ($X, $Y, $top) = @_; my %Z;
    for my $p (keys %$X) { for my $q (keys %$Y) { my $r = $p + $q; next if $r > $top; my $pr = mul($X->{$p}, $Y->{$q});
      $Z{$r} = defined $Z{$r} ? add($Z{$r}, $pr) : $pr } } \%Z };
  my %Ah = map { $_ => cst($A[$_]) } 0..$ZT;
  my $FE = $lmul->(\%F, \%E, $ZT);
  my %Psi = map { $_ => mul($FE->{$_}, $alg) } keys %$FE;            # F Exp[-dphi/z] alg
  my $okn = 1; for my $p (keys %Psi) { $okn = 0 if $p < 0 && !iszero($Psi{$p}) }
  chk($okn, "x=$x: [tau^k](F Exp[-dphi/z] alg) has no negative z-powers, k < $K  [Cor 7.1]");
  my $RHS = $lmul->(\%Psi, \%Ah, $M);
  my $okr = 1; for my $m (0..$M) { $okr = 0 unless iszero(sbt($RHS->{$m}, $rpt[$m])) }
  chk($okr, "x=$x: code Wick row (y = log U, tau-series) = h_w c/(c+z) F Exp[-dphi/z] alg  mod (tau^$K, z^" . ($M+1) . ")  [Prop 3.1 + Thm 5.2 + Prop 6.1]");
  # frame B: Exp[-ex] r^pt ; compare with (c/(c+z)) F Exp[-dphi/z] alg
  my @rB = map { my $m = $_; my $s = z0(); $s = add($s, scl($hwinv->[$_], $rpt[$m-$_])) for 0..$m; $s } 0..$M;
  my %cc = map { $_ => cst((R(-1)/$c)**$_) } 0..$ZT;
  my $RB = $lmul->(\%Psi, \%cc, $M);
  $okr = 1; for my $m (0..$M) { $okr = 0 unless iszero(sbt($RB->{$m}, $rB[$m])) }
  chk($okr, "x=$x: frame-B row Exp[-ex] r^pt = (c/(c+z)) F Exp[-dphi/z] alg  [Thm 7.2]");
  # convert to t: tau(t) = 1 - e^-t
  my $taut = z0(); $taut->[$_] = -((R(-1))**$_) / fact($_) for 1..$K-1;
  my @rBt = map { compose($_, $taut) } @rB;
  # rho from leaf_structure.md (data) and from tau/(t alg)
  my $rhoL = z0(); $rhoL->[0] = R(1); $rhoL->[1] = -(1-$x)**2/(2*$c**2); $rhoL->[2] = (1 - 8*$x + 15*$x**2 - 8*$x**3 + $x**4)/(6*$c**4);
  $rhoL->[3] = -(1-$x)**2 * (1 - 20*$x + 60*$x**2 - 20*$x**3 + $x**4)/(24*$c**6);
  $rhoL->[4] = (1 - 52*$x + 493*$x**2 - 1624*$x**3 + 2365*$x**4 - 1624*$x**5 + 493*$x**6 - 52*$x**7 + $x**8)/(120*$c**8);
  my $tovert = z0(); $tovert->[$_] = -((R(-1))**($_+1)) / fact($_+1) for 0..$K-1;   # tau/t = (1 - e^-t)/t
  my $rhoC = mul($tovert, inv(compose($alg, $taut)));
  my $okrho = 1; for my $i (0..4) { $okrho = 0 unless $rhoC->[$i] == $rhoL->[$i] }
  chk($okrho, "x=$x: rho of leaf_structure.md (t^0..t^4) = tau/(t alg)  [V1/Step 3 consistency with (1.4)]");
  # Delta_1 t^2/(-x c^3) = alg^2 t^2/tau^2 ; leaf_structure.md Step 0: 1 + (1-x)^2/c^2 t + (5-4x-6x^2-4x^3+5x^4)/(12c^4) t^2 + (1-x)^2/(12c^2) t^3
  my $d1n = mul(mul(compose($alg, $taut), compose($alg, $taut)), inv(mul($tovert, $tovert)));
  my @d1L = (R(1), (1-$x)**2/$c**2, (5 - 4*$x - 6*$x**2 - 4*$x**3 + 5*$x**4)/(12*$c**4), (1-$x)**2/(12*$c**2));
  my $okd = 1; for my $i (0..3) { $okd = 0 unless $d1n->[$i] == $d1L[$i] }
  chk($okd, "x=$x: -x delta/(tau^2 s*^2) = Delta_1 of leaf_structure.md Step 0 through t^1 (rel. t^3)");
  # dilaton: q_t(psi) = Sum_d (-1)^d r[d-1] psi^d, d = 1..M+1 ; rho q_t at t^k
  my @q = map { my $d = $_; scl(R(-1)**$d, $rBt[$d-1]) } 1..$M+1;     # q[d-1] = coefficient of psi^d (t-series)
  my @rq = map { mul($rhoL, $_) } @q;
  my $lam = sub { my ($j, $d) = @_; R($j)**$d / $c**($d+1) };          # [psi^d] 1/(c - j psi)
  my @delta = (undef,
    [ [1, -(1-$x)**2/2], [2, (1-$x)**2/2] ],
    [ [1, -(2 - 4*$x + 3*$x**2 - 4*$x**3 + 2*$x**4)/(12*$c**2)], [2, 6*(1-$x)**4/(12*$c**2)], [3, -(2-$x)**2*(1-2*$x)**2/(12*$c**2)] ],
    [ [1, -(1-$x)**2*(3 + 5*$x**2 + 3*$x**4)/(72*$c**4)], [2, (1-$x)**2*(7 - 20*$x + 22*$x**2 - 20*$x**3 + 7*$x**4)/(24*$c**4)],
      [3, -(1-$x)**2*(2-$x)**2*(1-2*$x)**2/(8*$c**4)], [4, (1-$x)**2*(3-$x)**2*(1-3*$x)**2/(36*$c**4)] ] );
  my $okq = 1;
  for my $d (1..$M+1) { my $q0 = ($d == 0) ? $c : -$c**2 * $lam->(1, $d); $okq = 0 unless $rq[$d-1][0] == $q0 }
  chk($okq, "x=$x: t^0: rho q_t = q_0 = c - c^2 lam_1 (psi^1..psi^" . ($M+1) . ")  [F2]");
  for my $k (1..3) { my $okk = 1;
    for my $d (1..$M+1) { my $v = R(0); $v += $_->[1] * $lam->($_->[0], $d) for @{$delta[$k]}; $okk = 0 unless $rq[$d-1][$k] == $v }
    chk($okk, "x=$x: t^$k: rho q_t = delta_$k of leaf_structure.md (psi^1..psi^" . ($M+1) . ")") }
  # R4: top coefficient of lam_{k+1} at order t^k: (-1)^(k+1) prod_{j=1}^k (k+1 - j c)^2 / (k! (k+1)! c^(2k-2))
  # extract it from the closed form: in z G = rho z r_{0,1}, residue at z = -c/(k+1) of F_k tau^k term only; equivalently,
  # the coefficient of lam_{k+1} is determined by the psi-series of [t^k] rho q_t as the dominant growth (k+1)^d / c^(d+1);
  # here we instead fit [t^k] rho q_t on lam_1..lam_{k+1} from psi^1..psi^(k+1) and test the remaining orders.
  for my $k (1..$M-1) { my $n = $k + 1;
    # solve Sum_{j=1}^{n} e_j lam_j[d] = rq[d-1][k], d = 1..n  (Vandermonde, exact Gauss elimination)
    my @Mx = map { my $d = $_; [ (map { $lam->($_, $d) } 1..$n), $rq[$d-1][$k] ] } 1..$n;
    for my $col (0..$n-1) { my ($piv) = grep { !$Mx[$_][$col]->is_zero } $col..$n-1; @Mx[$col, $piv] = @Mx[$piv, $col];
      for my $r (0..$n-1) { next if $r == $col; my $f = $Mx[$r][$col] / $Mx[$col][$col]; $Mx[$r][$_] -= $f * $Mx[$col][$_] for $col..$n } }
    my @e = map { $Mx[$_][$n] / $Mx[$_][$_] } 0..$n-1;
    my $okf = 1; for my $d ($n+1..$M+1) { my $v = R(0); $v += $e[$_-1] * $lam->($_, $d) for 1..$n; $okf = 0 unless $v == $rq[$d-1][$k] }
    my $pr = R(1); $pr *= ($k + 1 - $_ * $c)**2 for 1..$k;
    my $top = R(-1)**($k+1) * $pr / (fact($k) * fact($k+1) * $c**(2*$k-2));
    my $sum = R(0); $sum += $_ for @e;
    chk(($n < $M+1 ? $okf : 1) && $e[$n-1] == $top && $sum == 0, "x=$x: t^$k: [t^k] rho q_t in span(lam_1..lam_" . $n . ") (" . ($M+1-$n) . " extra orders), top coeff = (-1)^(k+1) prod (k+1-jc)^2/(k!(k+1)! c^(2k-2)), coefficients sum to 0") }
}
print "\n", ($bad ? "$bad FAILURES" : "all referee checks passed"), "\n";
