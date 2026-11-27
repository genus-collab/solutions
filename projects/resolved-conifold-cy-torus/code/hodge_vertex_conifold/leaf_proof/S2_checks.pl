#!/usr/bin/perl
# S2_checks.pl -- exact rational spot checks (Perl + Math::BigRat, no Mathematica) for section S2 of the leaf proof
# (leaf_proof/S2_metric_qde.md): Frobenius data at the torus point alpha0 = 0, alpha1 = eps1 = x, alpha2 = eps2 = 1,
# eigenvalues, metric Delta_i, rho = tau/(t alg), the conjugated PF identity (V7), the second QDE row, K_1 = 1.
# Usage: perl S2_checks.pl [N]      N = truncation order of the t- and tau-series (default 14).
# Part A: exact identities at rational points (x, tau) with rational delta = Sqrt[c^2 - 4 x tau].
# Part B: the operator identity (V7) on Laurent monomials tau^n, exactly, at rational z.
# Part C: the chain closed form (S1) -> H-row -> second QDE row -> K_1 = 1, as t-series at rational x, psi.
# Part D: the second row and the Wronskian form of K_1 = 1 from F = 2F1 alone, as tau-series at rational z.
# Part E: t-series of rho, Delta_1, dphi, Theta from the closed forms of S2 against the Mathematica outputs
#         (leaf_structure.md, hvc_leaf_proof_check_out.txt), at rational x.
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $N = @ARGV ? $ARGV[0] : 14;
my $M = $N + 6;                                   # internal length of all truncated series (indices 0..$M)
my ($ntest, $nfail) = (0, 0);
my %count;
sub ok { my ($cond, $msg) = @_; $ntest++; $count{(split /:/, $msg)[0]}++; unless ($cond) { $nfail++; print "FAIL: $msg\n"; } }

# ---------------------------------------------------------------- truncated power series (array refs 0..$M)
sub pz { [ map { R(0) } 0..$M ] }
sub pc { my $c = shift; my $a = pz(); $a->[0] = R($c); $a }
sub padd { my ($a, $b) = @_; [ map { $a->[$_] + $b->[$_] } 0..$M ] }
sub psub { my ($a, $b) = @_; [ map { $a->[$_] - $b->[$_] } 0..$M ] }
sub pscal { my ($s, $a) = @_; [ map { $a->[$_] * $s } 0..$M ] }
sub pmul { my ($a, $b) = @_; my $c = pz();
  for my $i (0..$M) { next if $a->[$i]->is_zero; for my $j (0..$M-$i) { $c->[$i+$j] += $a->[$i] * $b->[$j]; } } $c }
sub pinv { my $a = shift; die "pinv: zero constant term" if $a->[0]->is_zero; my $b = pz(); $b->[0] = 1 / $a->[0];
  for my $n (1..$M) { my $s = R(0); $s += $a->[$_] * $b->[$n-$_] for 1..$n; $b->[$n] = -$s / $a->[0]; } $b }
sub psqrt1 { my $a = shift; die "psqrt1: constant term != 1" unless $a->[0] == 1; my $b = pz(); $b->[0] = R(1);
  for my $n (1..$M) { my $s = R(0); $s += $b->[$_] * $b->[$n-$_] for 1..$n-1; $b->[$n] = ($a->[$n] - $s) / 2; } $b }
sub pexp0 { my $a = shift; die "pexp0: constant term != 0" unless $a->[0]->is_zero; my $e = pz(); $e->[0] = R(1);
  for my $n (1..$M) { my $s = R(0); $s += $a->[$_] * $e->[$n-$_] * $_ for 1..$n; $e->[$n] = $s / $n; } $e }
sub pder { my $a = shift; my $b = pz(); $b->[$_-1] = $a->[$_] * $_ for 1..$M; $b }      # top coefficient lost
sub pint { my $a = shift; my $b = pz(); $b->[$_] = $a->[$_-1] / $_ for 1..$M; $b }
sub pup { my $a = shift; my $b = pz(); $b->[$_] = $a->[$_-1] for 1..$M; $b }              # multiply by the variable
sub pcompose { my ($f, $g) = @_; die "pcompose: g(0) != 0" unless $g->[0]->is_zero;     # f(g)
  my $r = pz(); my $pw = pc(1); for my $k (0..$M) { $r = padd($r, pscal($f->[$k], $pw)); $pw = pmul($pw, $g); } $r }
sub pequal { my ($a, $b, $upto) = @_; for my $i (0..$upto) { return 0 unless $a->[$i] == $b->[$i]; } 1 }
sub pfirstdiff { my ($a, $b, $upto) = @_; for my $i (0..$upto) { return $i unless $a->[$i] == $b->[$i]; } -1 }

# ---------------------------------------------------------------- Laurent polynomials in tau (hash power => coeff)
sub lclean { my $a = shift; for my $k (keys %$a) { delete $a->{$k} if $a->{$k}->is_zero; } $a }
sub ladd { my ($a, $b) = @_; my %c = map { $_ => $a->{$_}->copy } keys %$a;
  for my $k (keys %$b) { $c{$k} = exists $c{$k} ? $c{$k} + $b->{$k} : $b->{$k}->copy; } lclean(\%c) }
sub lmono { my ($a, $s, $k) = @_; my %c; $c{$_ + $k} = $a->{$_} * $s for keys %$a; lclean(\%c) }   # s tau^k a
sub ldt { my $a = shift; my %c; for my $k (keys %$a) { $c{$k-1} = $a->{$k} * $k if $k != 0; } lclean(\%c) }  # d/dtau
sub lequal { my ($a, $b) = @_; my $d = ladd($a, lmono($b, R(-1), 0)); return !%$d; }

# ================================================================= part A
my @pts = ([R(2), R("5/8"), R(2)], [R("3/5"), R("7/15"), R("6/5")], [R(-3), R("5/12"), R(-3)]);
for my $p (@pts) {
  my ($x, $tau, $del) = @$p; my $c = 1 + $x;
  ok($del*$del == $c*$c - 4*$x*$tau, "A: delta^2 = c^2 - 4 x tau (x=$x)");
  my $Q = 1/(1 - $tau);
  my $s1 = 2/($c + $del); my $s2 = 1/($x*$tau*$s1);
  for my $s ($s1, $s2) { ok($x*$tau*$s*$s - $c*$s + 1 == 0, "A: x tau s^2 - c s + 1 = 0 (x=$x)"); }
  my $U1 = 1 - $tau*$s1; my $U2 = 1 - $tau*$s2;
  for my $U ($U1, $U2) { ok($x*$Q*$U*$U + $Q*(1-$x)*$U - 1 == 0, "A: x Q U^2 + Q(1-x) U - 1 = 0 (x=$x)"); }
  ok($U2 == 1 - 1/($x*$s1), "A: U2 = 1 - 1/(x s*)");
  ok($U1 + $U2 == ($x-1)/$x && $U1*$U2 == -1/($Q*$x), "A: Vieta in U (code SU, PU)");
  # U W'(U) = 0 at the critical points (W = x Log[1 - Q U] - Log[U] + Log[1 - U])
  for my $U ($U1, $U2) { ok(-$x*$Q*$U/(1 - $Q*$U) - 1 - $U/(1 - $U) == 0, "A: U W'(U) = 0 at U_i"); }
  my $LofU = sub { my $U = shift; -$Q*($x*$U + 1)/(1 - $Q) };            # rc_lib LofU at alpha0 = 0
  my ($L1, $L2) = ($LofU->($U1), $LofU->($U2));
  for my $pr ([$U1, $L1], [$U2, $L2]) { my ($U, $L) = @$pr;
    ok($L == 1/(1-$U) && $L == ($x*$U + 1)/$tau, "A: L(U) = 1/(1-U) = (xU+1)/tau on the critical locus");
    ok(-$x*$Q*$U/(1 - $Q*$U) == $L, "A: partial_t W = L on the critical locus"); }
  ok($L1 == 1/($tau*$s1) && $L1 == $c/$tau - $x*$s1 && $L2 == $x*$s1 && $L2 == 1/($tau*$s2), "A: L1 = 1/(tau s*) = c/tau - x s*, L2 = x s*");
  ok($L1 + $L2 == $c/$tau && $L1*$L2 == $x/$tau, "A: L1 + L2 = c/tau, L1 L2 = x/tau");
  ok($L1 - $L2 == $del/$tau, "A: L1 - L2 = delta/tau");
  # Frobenius structure: third derivatives of Prep = PrepCl + Li_3(Q)/(eps1 eps2) at the torus point
  my %P3 = ("000" => $c/$x**3, "001" => 1/$x**2, "011" => R(0), "111" => $Q/((1-$Q)*$x));
  my $P = sub { my @i = sort { $a <=> $b } @_; $P3{join("", @i)} };
  my @eta = ([$P->(0,0,0), $P->(0,0,1)], [$P->(0,1,0), $P->(0,1,1)]);
  ok($eta[0][0] == $c/$x**3 && $eta[0][1] == 1/$x**2 && $eta[1][0] == 1/$x**2 && $eta[1][1] == 0, "A: eta = {{c/x^3, 1/x^2}, {1/x^2, 0}}");
  my $det = $eta[0][0]*$eta[1][1] - $eta[0][1]*$eta[1][0];
  my @ei = ([$eta[1][1]/$det, -$eta[0][1]/$det], [-$eta[1][0]/$det, $eta[0][0]/$det]);
  ok($ei[0][0] == 0 && $ei[0][1] == $x**2 && $ei[1][0] == $x**2 && $ei[1][1] == -$c*$x, "A: eta^-1 = {{0, x^2}, {x^2, -c x}}");
  my %cch; for my $i (0,1) { for my $j (0,1) { for my $k (0,1) {
    my $s = R(0); $s += $ei[$i][$_] * $P->($_, $j, $k) for 0,1; $cch{"$i$j$k"} = $s; } } }
  ok($cch{"011"} == -$x/$tau && $cch{"111"} == $c/$tau, "A: cch: H*H = -(x/tau) 1 + (c/tau) H");
  ok($cch{"001"} == 0 && $cch{"101"} == 1 && $cch{"000"} == 1 && $cch{"100"} == 0, "A: cch: 1 is the unit");
  ok($Q/(1-$Q) == -1/$tau, "A: Q/(1-Q) = -1/tau");
  my $mult = sub { my ($a, $b) = @_; my @r;
    for my $i (0,1) { my $s = R(0); for my $j (0,1) { for my $k (0,1) { $s += $cch{"$i$j$k"} * $a->[$j] * $b->[$k]; } } push @r, $s; } \@r };
  my $veq = sub { my ($a, $b) = @_; $a->[0] == $b->[0] && $a->[1] == $b->[1] };
  # idempotents = columns of Inverse[Psi], Psi = {{1, L1}, {1, L2}}
  my $dd = $L2 - $L1; my @Pi = ([$L2/$dd, -$L1/$dd], [R(-1)/$dd, R(1)/$dd]);
  my $e1 = [$Pi[0][0], $Pi[1][0]]; my $e2 = [$Pi[0][1], $Pi[1][1]];
  ok($veq->($e1, [-$L2/($L1-$L2), 1/($L1-$L2)]) && $veq->($e2, [-$L1/($L2-$L1), 1/($L2-$L1)]), "A: e_1 = (H - L2)/(L1 - L2), e_2 = (H - L1)/(L2 - L1)");
  ok($veq->($mult->($e1, $e1), $e1) && $veq->($mult->($e2, $e2), $e2) && $veq->($mult->($e1, $e2), [R(0), R(0)]), "A: e_i e_j = delta_ij e_i");
  ok($veq->($mult->([R(0), R(1)], $e1), [$L1*$e1->[0], $L1*$e1->[1]]) && $veq->($mult->([R(0), R(1)], $e2), [$L2*$e2->[0], $L2*$e2->[1]]), "A: H e_i = L_i e_i");
  ok($e1->[0] + $e2->[0] == 1 && $e1->[1] + $e2->[1] == 0, "A: e_1 + e_2 = 1");
  my $etaf = sub { my ($a, $b) = @_; my $s = R(0); for my $i (0,1) { for my $j (0,1) { $s += $a->[$i]*$eta[$i][$j]*$b->[$j]; } } $s };
  my @ec; for my $i (0,1) { for my $j (0,1) { my $s = R(0); for my $a (0,1) { for my $b (0,1) { $s += $Pi[$a][$i]*$eta[$a][$b]*$Pi[$b][$j]; } } $ec[$i][$j] = $s; } }
  ok($ec[0][1] == 0 && $ec[1][0] == 0, "A: etacan diagonal");
  my ($D1, $D2) = (1/$ec[0][0], 1/$ec[1][1]);
  ok($etaf->($e1, $e1) == $etaf->($e1, [R(1), R(0)]) && $etaf->($e2, $e2) == $etaf->($e2, [R(1), R(0)]), "A: eta(e_i, e_i) = eta(e_i, 1)");
  my $W2 = sub { my $U = shift; -$x*$Q*$U/(1 - $Q*$U)**2 - $U/(1 - $U)**2 };           # (U d/dU)^2 W
  my $phi2 = sub { my $s = shift; -1/$s**2 - $x/(1 - $s)**2 + $tau**2/(1 - $tau*$s)**2 };  # phi_tau''(s)
  for my $q ([$D1, $L1, $L2, $U1, $s1], [$D2, $L2, $L1, $U2, $s2]) { my ($D, $Li, $Lj, $U, $s) = @$q;
    ok($D == $x**3*($Li - $Lj)/($x - $c*$Lj), "A: Delta_i = x^3 (L_i - L_j)/(x - c L_j)");
    ok($D == $x*$Li*(2*$x - $c*$Li), "A: Delta_i = x L_i (2x - c L_i)");
    ok($D == $x**3*($Li - $Lj)**2/($Lj*($c*$Lj - 2*$x)), "A: Delta_i = x^3 (L_i - L_j)^2/(L_j (c L_j - 2x)) (code quadratic form)");
    ok($D == $x**2 * $W2->($U), "A: Delta_i = x^2 W2(U_i)   (V1)");
    ok($W2->($U) == (1 - $tau*$s)**2 * $phi2->($s) / $tau**2, "A: W2(U_i) = (1 - tau s_i)^2 phi''(s_i)/tau^2");
    ok(-$x*$Q*$U/(1-$Q*$U)**2 == $Li*($x - $Li)/$x && $U/(1-$U)**2 == $Li*($Li - 1), "A: the two terms of W2 in terms of L"); }
  ok($D1 == -$x*$del/($tau**2*$s1**2) && $D2 == $x**3*$del*$s1**2, "A: Delta_1 = -x delta/(tau^2 s*^2), Delta_2 = x^3 delta s*^2");
  my $alg2 = (1 - $tau*$s1)**2 * $phi2->($s1) / (-$c**3/$x);
  ok((-1/(1/$c)**2 - $x/(1 - 1/$c)**2) == -$c**3/$x, "A: phi_0''(1/c) = -c^3/x");
  ok($alg2 == $del/($c**3*$s1**2) && $D1 == -$x*$c**3*$alg2/$tau**2, "A: alg^2 = delta/(c^3 s*^2), Delta_1 = -x c^3 alg^2/tau^2");
  ok(1/$D1 + 1/$D2 == $c/$x**3, "A: 1/Delta_1 + 1/Delta_2 = eta(1,1)");
  my $pair = sub { my ($a, $b) = @_; $x**2*($a->[0]*$b->[1] + $a->[1]*$b->[0]) - $c*$x*$a->[1]*$b->[1] };
  my ($f1, $f2) = ([R(1), $L1], [R(1), $L2]);
  ok($pair->($f1, $f2) == 0 && $pair->($f1, $f1) == $D1 && $pair->($f2, $f2) == $D2, "A: <(1,L_i),(1,L_j)> = delta_ij Delta_i");
  # A = {{0, 1}, {-x/tau, c/tau}} on covectors
  my @A = ([R(0), R(1)], [-$x/$tau, $c/$tau]);
  my $Av = sub { my $v = shift; [ $A[0][0]*$v->[0] + $A[0][1]*$v->[1], $A[1][0]*$v->[0] + $A[1][1]*$v->[1] ] };
  ok($veq->($Av->($f1), [$L1, $L1*$L1]) && $veq->($Av->($f2), [$L2, $L2*$L2]), "A: A (1, L_i) = L_i (1, L_i)");
  my @G = ([R(0), $x**2], [$x**2, -$c*$x]);
  my @GA; for my $i (0,1) { for my $j (0,1) { $GA[$i][$j] = $G[$i][0]*$A[0][$j] + $G[$i][1]*$A[1][$j]; } }
  ok($GA[0][1] == $GA[1][0], "A: eta^-1 A symmetric (A self-adjoint)");
  # also: A is the transpose of the matrix of H* in the flat basis
  my $H1 = $mult->([R(0), R(1)], [R(1), R(0)]); my $HH = $mult->([R(0), R(1)], [R(0), R(1)]);
  ok($A[0][0] == $H1->[0] && $A[0][1] == $H1->[1] && $A[1][0] == $HH->[0] && $A[1][1] == $HH->[1], "A: A = transpose of H* (rows = H*1, H*H)");
  # spectral decomposition A - L1 = (L2 - L1) P2, P2 a = pi(a)/Delta_2 (1, L2), pi(a) = <a, (1, L2)>
  for my $a ([R(1), R(0)], [R(0), R(1)], [R("2/3"), R("-5/7")]) {
    my $lhs = $Av->($a); $lhs = [$lhs->[0] - $L1*$a->[0], $lhs->[1] - $L1*$a->[1]];
    my $pi = $pair->($a, $f2); my $rhs = [($L2 - $L1)*$pi/$D2, ($L2 - $L1)*$pi*$L2/$D2];
    ok($veq->($lhs, $rhs), "A: A - L1 = (L2 - L1) P2");
    ok($pi == $x**2*$L2*$a->[0] + ($x**2 - $c*$x*$L2)*$a->[1], "A: pi(a) = x^2 L2 a0 + (x^2 - c x L2) a1");
    ok($pi == $x**3*$s1*($a->[0] - $tau*$s1*$a->[1]), "A: pi(a) = x^3 s* (a0 - tau s* a1)");
    for my $b ([R(1), R(0)], [R(0), R(1)], [R("-4/9"), R("3/2")]) {
      ok($pair->($lhs, $b) == -($x**3/$tau)*($a->[0] - $tau*$s1*$a->[1])*($b->[0] - $tau*$s1*$b->[1]), "A: <(A - L1) a, b> = -(x^3/tau)(a0 - tau s* a1)(b0 - tau s* b1)"); } }
  ok($x - $c*$L2 == -$x**2*$tau*$s1**2 && $x**2 - $c*$x*$L2 == -$x**3*$tau*$s1**2, "A: x - c L2 = -x^2 tau s*^2");
  my $Theta = ($L2 - $L1)*($c*$x*$L2 - $x**2)**2/$D2;
  ok($Theta == -$x**3*$tau*$s1**2 && $Theta == $x*$c*(1 - $x*$s1) - $x, "A: Theta = -x^3 tau s*^2 = x c dphi' - x");
  ok((1 - $tau)*$s1/(1 - $tau*$s1) == 1 - $x*$s1 && $L1 - (1 - $x*$s1) == $c/$tau - 1, "A: dphi' = (1-tau)s*/(1-tau s*) = 1 - x s*, L1 - dphi' = c/tau - 1");
}

# ================================================================= part B: (V7) in conjugated form, exactly
# (D + e)^2 f - (c/tau)(D + e) f + (x/tau) f = ((1 - tau) z^2/tau) Gauss[f],  D = z (1 - tau) d/dtau, e = (c+z)/tau - (1+z)
for my $x (R(2), R("3/5"), R(-3), R("-7/4")) { my $c = 1 + $x;
  for my $z (R("2/7"), R("-5/3"), R("4/9"), R(11)) {
    my ($al, $ga) = (1/$z + 1, $c/$z + 2);
    my $D = sub { my $f = shift; my $g = ldt($f); ladd(lmono($g, $z, 0), lmono($g, -$z, 1)) };
    my $De = sub { my $f = shift; ladd($D->($f), ladd(lmono($f, $c + $z, -1), lmono($f, -(1 + $z), 0))) };
    my $gauss = sub { my $f = shift; my $f1 = ldt($f); my $f2 = ldt($f1);
      my $r = ladd(lmono($f2, R(1), 1), lmono($f2, R(-1), 2)); $r = ladd($r, lmono($f1, $ga, 0));
      $r = ladd($r, lmono($f1, -(2*$al + 1), 1)); ladd($r, lmono($f, -$al*$al, 0)) };
    my @tests = map { { $_ => R(1) } } (-3..7);
    push @tests, { -2 => R("3/5"), 0 => R(-2), 1 => R("7/3"), 4 => R("1/11") };
    for my $f (@tests) {
      my $g1 = $De->($f); my $lhs = ladd(ladd($De->($g1), lmono($g1, -$c, -1)), lmono($f, $x, -1));
      my $G = $gauss->($f); my $rhs = ladd(lmono($G, $z*$z, -1), lmono($G, -$z*$z, 0));
      ok(lequal($lhs, $rhs), "B: (V7) conjugated identity on tau^n (x=$x, z=$z)"); }
    # D(P)/P = e for P = (1-tau)^(-x/z) tau^(c/z+1): z (1-tau) [ (x/z)/(1-tau) + (c/z+1)/tau ] = x + (c+z)(1-tau)/tau
    ok($x + ($c + $z)*(1 - R("1/3"))/R("1/3") == ($c + $z)/R("1/3") - (1 + $z), "B: e = D(P)/P at tau = 1/3");
  } }

# ================================================================= helpers for parts C, D, E
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
sub Fcoef { my ($z, $c) = @_; my @F;                   # [tau^n] 2F1(1/z+1, 1/z+1; c/z+2; tau), product formula of (V3)
  for my $n (0..$M) { my $num = R(1); my $den = fact($n) * $z**$n;
    $num *= (1 + $_*$z)**2 for 1..$n; $den *= ($c + $_*$z) for 2..$n+1; push @F, $num/$den; } \@F }
sub goodz { my ($c, @cand) = @_; grep { my $r = $c/$_; !($r->is_int && abs($r) <= $M + 3) } @cand }

# ================================================================= part D: second row and K_1 from F alone (tau-series)
# Psi1 := tau (D + e) F = z tau (1-tau) F' + (c + z - (1+z) tau) F;  row 2  <=>  z (1-tau) Psi1' = Psi1 - x F;
# K_1 = 1  <=>  tau^2 <Phi(z), Phi(-z)> = -x c (c^2 - z^2), Phi = (F, Psi1/tau).
for my $x (R(2), R("3/5"), R(-3), R("-7/4")) { my $c = 1 + $x;
  for my $z (goodz($c, R("2/7"), R("-5/3"), R("4/9"), R("3/7"), R("5/11"))) {
    my %S; for my $zz ($z, -$z) {
      my $F = Fcoef($zz, $c); my ($al, $ga) = (1/$zz + 1, $c/$zz + 2);
      for my $n (0..$M-1) { ok(($n+1)*($n+$ga)*$F->[$n+1] == ($n+$al)**2*$F->[$n], "D: Gauss recursion of the product formula"); last if $n > 8; }
      my $Fp = pder($F); my $tau = pup(pc(1)); my $omt = psub(pc(1), $tau);
      my $Psi1 = padd(pscal($zz, pmul(pmul($tau, $omt), $Fp)), pmul(padd(pc($c + $zz), pscal(-(1 + $zz), $tau)), $F));
      my $row2 = psub(pscal($zz, pmul($omt, pder($Psi1))), psub($Psi1, pscal($x, $F)));
      ok(pequal($row2, pz(), $N), "D: second QDE row z(1-tau) Psi1' = Psi1 - x F to tau^$N (x=$x, z=$zz)");
      $S{$zz} = [$F, $Psi1]; }
    my ($F, $P1) = @{$S{$z}}; my ($Fm, $P1m) = @{$S{-$z}};
    my $W = psub(pscal($x**2, pup(padd(pmul($F, $P1m), pmul($P1, $Fm)))), pscal($c*$x, pmul($P1, $P1m)));
    ok(pequal($W, pc(-$x*$c*($c**2 - $z**2)), $N), "D: tau^2 <Phi(z), Phi(-z)> = -x c (c^2 - z^2) to tau^$N (x=$x, z=$z)");
  } }

# ================================================================= parts C, E: t-series at rational x
for my $x (R(2), R("3/5"), R(-3)) { my $c = 1 + $x;
  my $tt = pz(); $tt->[$_] = R((($_ % 2) ? -1 : 1)) / fact($_ + 1) for 0..$M;      # tau/t = (1 - e^-t)/t
  my $tau = pup($tt); my $toverTau = pinv($tt);
  my $dc = psqrt1(psub(pc(1), pscal(4*$x/$c**2, $tau)));                             # delta/c
  my $cs = pscal(R(2), pinv(padd(pc(1), $dc)));                                      # c s*
  my $ss = pscal(1/$c, $cs);                                                         # s*
  ok(pequal(padd(psub(pscal($x, pmul($tau, pmul($ss, $ss))), pscal($c, $ss)), pc(1)), pz(), $N), "C: x tau s*^2 - c s* + 1 = 0 (series, x=$x)");
  my $dphip = psub(pc(1), pscal($x, $ss)); my $dphi = pint($dphip);
  my $omt = psub(pc(1), $tau);
  ok(pequal(pmul(pmul($omt, $ss), pinv(psub(pc(1), pmul($tau, $ss)))), $dphip, $N), "C: dphi' = (1-tau) s*/(1 - tau s*) = 1 - x s* (series, x=$x)");
  my $alg = pmul(psqrt1($dc), pinv($cs));
  my $rho = pmul($tt, pinv($alg));
  my $D1n = pmul(pmul($toverTau, $toverTau), pmul($alg, $alg));                      # Delta_1 t^2/(-x c^3)
  ok(pequal($D1n, pinv(pmul($rho, $rho)), $N), "C: rho^-2 = Delta_1 t^2/(-x c^3) (x=$x)");
  my $D2 = pscal($x**3*$c, pmul($dc, pmul($ss, $ss)));                              # Delta_2 = x^3 delta s*^2
  ok(pequal(psub($D2, psub(pc($x**3/$c), pscal($x**5/$c**5, pmul($tau, $tau)))), pz(), 2), "C: Delta_2 = x^3/c - x^5 tau^2/c^5 + O(tau^3) (x=$x)");
  my $Theta = pscal(-$x**3, pmul($tau, pmul($ss, $ss)));
  ok(pequal($Theta, psub(pscal($x*$c, $dphip), pc($x)), $N), "C: Theta = x c dphi' - x (series, x=$x)");
  # ---- part E: comparison with the Mathematica outputs
  my $one = R(1);
  my @rhoM = ($one, -(1-$x)**2/(2*$c**2), (1 - 8*$x + 15*$x**2 - 8*$x**3 + $x**4)/(6*$c**4),
     -(1-$x)**2*(1 - 20*$x + 60*$x**2 - 20*$x**3 + $x**4)/(24*$c**6),
     (1 - 52*$x + 493*$x**2 - 1624*$x**3 + 2365*$x**4 - 1624*$x**5 + 493*$x**6 - 52*$x**7 + $x**8)/(120*$c**8));
  ok(pequal($rho, [@rhoM, map { R(0) } 5..$M], 4), "E: rho = tau/(t alg) matches leaf_structure.md to t^4 (x=$x)");
  my @D1M = ($one, (1-$x)**2/$c**2, (5 - 4*$x - 6*$x**2 - 4*$x**3 + 5*$x**4)/(12*$c**4), (1-$x)**2/(12*$c**2));
  ok(pequal($D1n, [@D1M, map { R(0) } 4..$M], 3), "E: Delta_1 t^2/(-x c^3) matches leaf_structure.md to t^3 (x=$x)");
  my $p4 = 1 - 8*$x + 12*$x**2 - 8*$x**3 + $x**4; my $p5 = 1 - 20*$x + 42*$x**2 - 20*$x**3 + $x**4;
  my $p8 = 1 - 52*$x + 418*$x**2 - 1204*$x**3 + 1690*$x**4 - 1204*$x**5 + 418*$x**6 - 52*$x**7 + $x**8;
  my $q8 = 1 - 112*$x + 1528*$x**2 - 5824*$x**3 + 8830*$x**4 - 5824*$x**5 + 1528*$x**6 - 112*$x**7 + $x**8;
  my $q12 = 1 - 240*$x + 6576*$x**2 - 56480*$x**3 + 223095*$x**4 - 483888*$x**5 + 621600*$x**6 - 483888*$x**7 + 223095*$x**8 - 56480*$x**9 + 6576*$x**10 - 240*$x**11 + $x**12;
  my @dM = (R(0), 1/$c, -$x**2/(2*$c**3), ($x-1)**2*$x**2/(6*$c**5), -$x**2*$p4/(24*$c**7), ($x-1)**2*$x**2*$p5/(120*$c**9),
     -$x**2*$p8/(720*$c**11), ($x-1)**2*$x**2*$q8/(5040*$c**13));
  ok(pequal($dphi, [@dM, map { R(0) } 8..$M], 7), "E: dphi = Int (1 - x s*) dt matches hvc_leaf_proof_check_out.txt to t^7 (x=$x)");
  my @thM = (R(0), -$x**3/$c**2, ($x-1)**2*$x**3/(2*$c**4), -$x**3*$p4/(6*$c**6), ($x-1)**2*$x**3*$p5/(24*$c**8),
     -$x**3*$p8/(120*$c**10), ($x-1)**2*$x**3*$q8/(720*$c**12), -$x**3*$q12/(5040*$c**14));
  ok(pequal($Theta, [@thM, map { R(0) } 8..$M], 7), "E: Theta = -x^3 tau s*^2 matches hvc_leaf_proof_check_out.txt to t^7 (x=$x)");
  # ---- part C: closed form -> H-row -> second row -> K_1 = 1
  my $tL1 = psub(pscal($c, $toverTau), pscal($x, pup($ss)));                         # t L_1 = c t/tau - x t s*
  my $L2 = pscal($x, $ss);
  my $Gof = sub { my $z = shift;                                                     # G(z) = c/(c+z) (tau/t) F(z, tau) Exp[-dphi/z]
    my $Ft = pcompose(Fcoef($z, $c), $tau); my $ex = pexp0(pscal(-1/$z, $dphi));
    pscal($c/($c + $z), pmul($tt, pmul($Ft, $ex))) };
  for my $psi (goodz($c, R("2/7"), R("-5/3"), R("4/9"), R("3/7"))) {
    my %g; my %Sig;
    for my $pp ($psi, -$psi) {
      my $g = $Gof->(-$pp);                                                          # g(pp) = G(-pp)
      my $Sig = psub(pmul(psub($tL1, pc($pp)), $g), pscal($pp, pup(pder($g))));     # t sigma_1(pp) = (t L1 - pp) g - pp t g'
      ok($g->[0] == $c/($c - $pp), "C: F2: g(psi; 0) = c/(c - psi) (x=$x, psi=$pp)");
      ok($Sig->[0] == $c, "C: [t^-1] sigma_1 = c (x=$x, psi=$pp)");
      my $row2 = psub(padd(pscal($pp, pder($Sig)), pmul($L2, $Sig)), pscal($x, pmul($toverTau, $g)));
      ok(pequal($row2, pz(), $N), "C: second row for sigma: psi (t sigma_1)' + L2 (t sigma_1) - x (t/tau) g = 0 to t^$N (x=$x, psi=$pp)");
      # t pi(sigma) = x^2 L2 t g + (x^2 - c x L2) (t sigma_1)  versus  -(x^2 - c x L2) psi (g + t g')
      my $m2 = psub(pc($x**2), pscal($c*$x, $L2));
      my $tpi = padd(pscal($x**2, pmul($L2, pup($g))), pmul($m2, $Sig));
      ok(pequal($tpi, pscal(-$pp, pmul($m2, padd($g, pup(pder($g))))), $N), "C: pi(sigma) = -(x^2 - c x L2)(psi/t) d/dt(t g) to t^$N (x=$x, psi=$pp)");
      $g{$pp} = $g; $Sig{$pp} = $Sig; }
    my $K = psub(pscal($x**2, pup(padd(pmul($g{$psi}, $Sig{-$psi}), pmul($Sig{$psi}, $g{-$psi})))), pscal($c*$x, pmul($Sig{$psi}, $Sig{-$psi})));
    ok(pequal($K, pc(-$x*$c**3), $N), "C: t^2 <sigma(psi), sigma(-psi)> = -x c^3 (K_1 = 1) to t^$N (x=$x, psi=$psi)");
  } }

print "S2_checks.pl, N = $N: $ntest checks, $nfail failures\n";
for my $k (sort keys %count) { print "  $k: $count{$k}\n"; }
