#!/usr/bin/perl
# ref_indep.pl -- referee's independent recomputation (Perl + Math::BigRat) of the frame-B factors of rc_lib.m / hvc_lib.m
# at a rational x, directly from the code's DEFINITIONS:
#   * point-frame identity rows = wickDirect polynomials in W_k = (U d/dU)^k spotU at U1(t), U2(t)  (no closed form used);
#   * Delta_i = 1/etacan[[i,i]] with Psi = {{1,L1},{1,L2}}, L = LofU;
#   * H-row of rmatrixCore; twistRows[hCoeffs[{c,-x,-1}, -1]];
#   * EdN of evalQ (formula (1.1) of S4).
# Then it tests the claims of leaf_proof/S4_legs_closed.md and compares with the Mathematica leaf tables.
use strict; use warnings; use Math::BigRat;

our $NT = 8;            # t-series kept to t^NT (exactness of derived quantities is tracked below)
my $MZ = 6;             # rows r[0..MZ]
my @XS = @ARGV ? @ARGV : ("3/7");
use File::Basename (); use File::Spec ();
my $ROOT = File::Spec->rel2abs(File::Basename::dirname(__FILE__)) . "/../../..";

sub R { my $v = shift; return ref($v) ? $v->copy : Math::BigRat->new("$v"); }
sub sz { return [ map { Math::BigRat->new(0) } 0..$NT ]; }
sub sc { my $c = shift; my $s = sz(); $s->[0] = R($c); return $s; }
sub sadd { my ($a, $b) = @_; return [ map { $a->[$_] + $b->[$_] } 0..$NT ]; }
sub ssub { my ($a, $b) = @_; return [ map { $a->[$_] - $b->[$_] } 0..$NT ]; }
sub ssc  { my ($k, $a) = @_; my $kk = R($k); return [ map { $a->[$_] * $kk } 0..$NT ]; }
sub smul { my ($a, $b) = @_; my $r = sz();
  for my $i (0..$NT) { next if $a->[$i]->is_zero; for my $j (0..$NT-$i) { next if $b->[$j]->is_zero; $r->[$i+$j] += $a->[$i] * $b->[$j]; } }
  return $r; }
sub sinv { my $a = shift; die "sinv: zero constant term" if $a->[0]->is_zero; my $r = sz(); $r->[0] = R(1) / $a->[0];
  for my $n (1..$NT) { my $s = R(0); for my $k (1..$n) { $s += $a->[$k] * $r->[$n-$k]; } $r->[$n] = -$s / $a->[0]; } return $r; }
sub sdiv { my ($a, $b) = @_; return smul($a, sinv($b)); }
sub sder { my $a = shift; my $r = sz(); for my $n (1..$NT) { $r->[$n-1] = $a->[$n] * $n; } return $r; }
sub sint { my $a = shift; my $r = sz(); for my $n (0..$NT-1) { $r->[$n+1] = $a->[$n] / ($n+1); } return $r; }
sub sexp { my $a = shift; die "sexp: nonzero constant" unless $a->[0]->is_zero; my $r = sz(); $r->[0] = R(1);
  for my $n (1..$NT) { my $s = R(0); for my $k (1..$n) { $s += $a->[$k] * $r->[$n-$k] * $k; } $r->[$n] = $s / $n; } return $r; }
sub slog { my $a = shift; die "slog: constant term not 1" unless $a->[0] == 1; return sint(smul(sder($a), sinv($a))); }
sub spow { my ($a, $p) = @_; my $a0 = $a->[0]; die "spow: constant term not 1 ($a0)" unless $a0 == 1; return sexp(ssc($p, slog($a))); }
sub sup  { my ($a, $m) = @_; my $r = sz(); for my $n ($m..$NT) { $r->[$n] = $a->[$n-$m]->copy; } return $r; }
sub sdown { my ($a, $m, $lab) = @_; for my $n (0..$m-1) { die "sdown($lab): nonzero coefficient t^$n = $a->[$n]" unless $a->[$n]->is_zero; }
  my $r = sz(); for my $n (0..$NT-$m) { $r->[$n] = $a->[$n+$m]->copy; } return $r; }
sub scompose { my ($f, $g) = @_; die "compose" unless $g->[0]->is_zero; my $r = sz(); my $p = sc(1);
  for my $n (0..$NT) { $r = sadd($r, ssc($f->[$n], $p)); $p = smul($p, $g); } return $r; }
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 1..$n; return $f; }
sub hyp21 { my ($a, $b, $g) = @_; my @f = (R(1)); for my $n (1..$NT) { push @f, $f[-1] * ($a + $n - 1) * ($b + $n - 1) / (($g + $n - 1) * $n); } return \@f; }

# partitions of N into exactly P parts, each >= min (nondecreasing lists)
sub partitions { my ($N, $P, $min) = @_;
  return ([]) if $P == 0 && $N == 0;
  return () if $P == 0 || $N < $P * $min;
  my @out; for my $first ($min .. int($N / $P)) { for my $rest (partitions($N - $first, $P - 1, $first)) { push @out, [$first, @$rest]; } }
  return @out; }
# wickZ[[m+1]] as a list of [coef, n, {k => mult}] : coef * W2^-n * prod W_k^mult   (wickDirect, rc_lib.m s.3)
sub wickTerms { my $m = shift; return ([R(1), 0, {}]) if $m == 0; my @T;
  for my $n ($m + 1 .. 3 * $m) { my $j = $n - $m;
    for my $p (partitions(2 * $n, $j, 3)) { my %mu; $mu{$_}++ for @$p;
      my $co = fact(2 * $n) * R((-1) ** $n) / (fact($n) * R(2) ** $n);
      for my $k (keys %mu) { $co /= (fact($k) ** $mu{$k}) * fact($mu{$k}); }
      push @T, [$co, $n, \%mu]; } }
  return @T; }
# P_k(Y) := theta^k log(1 - a U) as a polynomial in Y = (1 - a U)^-1:  P_1 = 1 - Y, theta Y^n = n (Y^(n+1) - Y^n)
my @Ppoly; $Ppoly[1] = [R(1), R(-1)];
for my $k (2..2 * $MZ + 2) { my $p = $Ppoly[$k-1]; my @q = map { R(0) } 0..$k;
  for my $n (1..$#$p) { next if $p->[$n]->is_zero; $q[$n+1] += $p->[$n] * $n; $q[$n] -= $p->[$n] * $n; } $Ppoly[$k] = \@q; }

my (%nC, %nB);
sub chk { my ($cat, $lab, $a, $b) = @_; $nC{$cat}++; if ($a != $b) { $nB{$cat}++; print "  MISMATCH [$cat] $lab: $a vs $b\n" if ($nB{$cat} // 0) <= 6; } }

# Mathematica tables (InputForm), evaluated at rational x
sub splitTop { my $s = shift; my @out; my $d = 0; my $cur = "";
  for my $ch (split //, $s) { if ($ch eq '(' || $ch eq '{' || $ch eq '[') { $d++; } elsif ($ch eq ')' || $ch eq '}' || $ch eq ']') { $d--; }
    if ($ch eq ',' && $d == 0) { push @out, $cur; $cur = ""; } else { $cur .= $ch; } }
  push @out, $cur if $cur =~ /\S/; s/^\s+|\s+$//g for @out; return @out; }
sub evalMma { my ($expr, $xv) = @_; my $p = $expr; $p =~ s/\^/**/g; $p =~ s/\bx\b/\$x/g;
  my $x = $xv; my $v; { use bigrat; $v = eval $p; } die "eval failed ($@): $p" if $@ || !defined $v; return Math::BigRat->new("$v"); }
sub slurp { my $f = shift; open(my $fh, "<", $f) or die "cannot open $f"; local $/; my $s = <$fh>; close $fh; return $s; }
my %leg;
for my $f ("$ROOT/hvc_leaf_structure_out.txt", "$ROOT/hvc_leaf_structure_out_kt5.txt") { for my $line (split /\n/, slurp($f)) {
  if ($line =~ /^\(ii\) leg j = (\d+), b = (\d+), t\^(\d+): lam_0\.\.lam_\d+ = \{(.*)\}; checks/) { $leg{"$1,$2,$3"} = [ splitTop($4) ]; } } }

my @B2 = (undef, R("1/6"), R("-1/30"), R("1/42"), R("-1/30"));
my $t0 = time;
for my $xs (@XS) {
  my $x = R($xs); my $c = 1 + $x; print "==================== x = $x (c = $c)\n";
  # ---- t-series: tau, Q, delta, s*, U1, U2
  my $tot = sz(); $tot->[$_] = R((-1) ** $_) / fact($_ + 1) for 0..$NT;     # tau/t
  my $tau = sup($tot, 1);
  my $Qs = sz(); $Qs->[$_] = R(1) / fact($_) for 0..$NT;                       # e^t
  my $em1t = sz(); $em1t->[$_] = R(1) / fact($_ + 1) for 0..$NT;               # (e^t - 1)/t
  my $tem1 = sinv($em1t);                                                       # t/(e^t - 1)
  my $delta = ssc($c, spow(ssub(sc(1), ssc(4 * $x / ($c * $c), $tau)), R("1/2")));
  my $sst = ssc(2, sinv(sadd(sc($c), $delta)));
  my $U1 = ssub(sc(1), smul($tau, $sst));
  my $U2 = ssub(sc(1), ssc(R(1) / (2 * $x), sadd(sc($c), $delta)));
  # sanity: both roots of Q x U^2 + Q (1-x) U - 1
  for my $U ($U1, $U2) { my $e = ssub(smul($Qs, sadd(ssc($x, smul($U, $U)), ssc(1 - $x, $U))), sc(1)); chk("roots", "Leq", $e->[$_], R(0)) for 0..$NT; }
  # ---- W_k at the two roots. Root 1 normalised: w_k = t^k W_k (a = t/(1-U1), b = t/(1-Q U1)); root 2 plain.
  my $a1 = sdiv(sinv($tot), $sst);                                              # t/(tau s*)
  my $b1 = sdiv($tem1, ssub($sst, sc(1)));                                      # t/((e^t-1)(s*-1)) ; 1 - Q U1 = (e^t-1)(s*-1)
  my $A2 = sinv(ssub(sc(1), $U2)); my $B2 = sinv(ssub(sc(1), smul($Qs, $U2)));
  my (%W1n, %W2);
  for my $k (1..2 * $MZ + 2) { my $p = $Ppoly[$k];
    my $s1 = sz(); my $s2 = sz(); my ($pa, $pb, $qa, $qb) = (sc(1), sc(1), sc(1), sc(1));
    for my $n (0..$#$p) { if ($n > 0) { $pa = smul($pa, $a1); $pb = smul($pb, $b1); $qa = smul($qa, $A2); $qb = smul($qb, $B2); }
      next if $p->[$n]->is_zero;
      $s1 = sadd($s1, ssc($p->[$n], sup(sadd(ssc($x, $pb), $pa), $k - $n)));
      $s2 = sadd($s2, ssc($p->[$n], sadd(ssc($x, $qb), $qa))); }
    if ($k == 1) { $s1 = ssub($s1, sup(sc(1), 1)); $s2 = ssub($s2, sc(1)); }
    $W1n{$k} = $s1; $W2{$k} = $s2; }
  chk("W1 = 0 at roots", "root1", $W1n{1}->[$_], R(0)) for 0..$NT; chk("W1 = 0 at roots", "root2", $W2{1}->[$_], R(0)) for 0..$NT;
  # ---- Wick rows (point frame)
  my (@r01pt, @r02pt);
  for my $root (1, 2) { my $W = $root == 1 ? \%W1n : \%W2; my $iw2 = sinv($W->{2}); my (%pw, @ipw); $ipw[0] = sc(1); $ipw[$_] = smul($ipw[$_-1], $iw2) for 1..3 * $MZ;
    my $pow; $pow = sub { my ($k, $e) = @_; my $key = "$k,$e"; return $pw{$key} if exists $pw{$key};
      my $s = ($e == 0) ? sc(1) : smul($W->{$k}, $pow->($k, $e - 1)); $pw{$key} = $s; return $s; };
    for my $m (0..$MZ) { my $s = sz();
      for my $T (wickTerms($m)) { my ($co, $n, $mu) = @$T; my $v = ssc($co, $ipw[$n]); for my $k (keys %$mu) { $v = smul($v, $pow->($k, $mu->{$k})); } $s = sadd($s, $v); }
      if ($root == 1) { $r01pt[$m] = $s; } else { $r02pt[$m] = $s; } } }
  printf "  Wick rows done (%d s)\n", time - $t0;
  # ---- L_i (LofU), Delta_i = 1/etacan[[i,i]] = x^3 (L_j - L_i)^2/(L_j (c L_j - 2x))
  my $tL1 = smul(smul($Qs, sadd(ssc($x, $U1), sc(1))), $tem1);                 # t L1 = Q (x U1 + 1) t/(Q - 1)   (exact to NT)
  my $tL2 = smul(smul($Qs, sadd(ssc($x, $U2), sc(1))), $tem1);
  my $L2 = sdown($tL2, 1, "L2");                                               # exact to NT-1
  my $t2D1 = ssc($x ** 3, sdiv(smul(ssub($tL2, $tL1), ssub($tL2, $tL1)), smul($L2, ssub(ssc($c, $L2), sc(2 * $x)))));
  my $D2 = ssc($x ** 3, sdiv(smul(ssub($tL1, $tL2), ssub($tL1, $tL2)), smul($tL1, ssub(ssc($c, $tL1), sup(sc(2 * $x), 1)))));
  chk("Delta_1 leading", "t^2 Delta_1 (0) = -x c^3", $t2D1->[0], -$x * $c ** 3);
  chk("Delta_2 (10.1)", "Delta_2 = x^3 delta s*^2, t^$_", $D2->[$_], ssc($x ** 3, smul($delta, smul($sst, $sst)))->[$_]) for 0..$NT;
  chk("Delta_2 (10.1)", "Delta_2(0) = x^3/c", $D2->[0], $x ** 3 / $c); chk("Delta_2 (10.1)", "Delta_2'(0) = 0", $D2->[1], R(0));
  my $dl1 = ssc(R("1/2"), sdiv(sder($t2D1), $t2D1));                          # Delta_1'/(2 Delta_1) = dl1 - 1/t
  my $dl2 = ssc(R("1/2"), sdiv(sder($D2), $D2));
  my $rho = spow(ssc(R(-1) / ($x * $c ** 3), $t2D1), R("-1/2"));
  # V1 check: Delta_i = x^2 W_2(U_i)
  chk("V1", "t^2 Delta_1 = x^2 w_2, t^$_", $t2D1->[$_], $x * $x * $W1n{2}->[$_]) for 0..$NT - 1;
  chk("V1", "Delta_2 = x^2 W_2(U2), t^$_", $D2->[$_], $x * $x * $W2{2}->[$_]) for 0..$NT;
  # ---- H-rows (point frame): r12 = L2 r02 + r02[m-1]' - dl2 r02[m-1];  T1 := t r11 = tL1 r01 + r01[m-1] + t (r01[m-1]' - dl1 r01[m-1])
  my (@r12pt, @T1pt);
  for my $m (0..$MZ) {
    $r12pt[$m] = smul($L2, $r02pt[$m]); $T1pt[$m] = smul($tL1, $r01pt[$m]);
    if ($m >= 1) { $r12pt[$m] = sadd($r12pt[$m], ssub(sder($r02pt[$m-1]), smul($dl2, $r02pt[$m-1])));
      $T1pt[$m] = sadd($T1pt[$m], sadd($r01pt[$m-1], sup(ssub(sder($r01pt[$m-1]), smul($dl1, $r01pt[$m-1])), 1))); } }
  # ---- frame B: multiply every row by h_w^-1 = exp(+Sum beta_k S_(2k-1) z^(2k-1)),  S_m = c^-m - x^-m - 1
  my @ex = map { R(0) } 0..$MZ;
  for my $k (1..4) { my $m = 2 * $k - 1; last if $m > $MZ; $ex[$m] = $B2[$k] / (2 * $k * (2 * $k - 1)) * (R(1) / $c ** $m - R(1) / $x ** $m - 1); }
  my $zexp = sub { my ($e) = @_; my @r = (R(1)); for my $n (1..$MZ) { my $s = R(0); for my $k (1..$n) { $s += $e->[$k] * $r[$n-$k] * $k; } push @r, $s / $n; } return \@r; };
  my $hinv = $zexp->(\@ex); my $h2 = $zexp->([ map { $_ * -2 } @ex ]); my $hm2 = $zexp->([ map { $_ * 2 } @ex ]);
  my $tw = sub { my ($rows) = @_; my @o; for my $m (0..$MZ) { my $s = sz(); $s = sadd($s, ssc($hinv->[$m - $_], $rows->[$_])) for 0..$m; $o[$m] = $s; } return \@o; };
  my $r01 = $tw->(\@r01pt); my $T1 = $tw->(\@T1pt); my $r02 = $tw->(\@r02pt); my $r12 = $tw->(\@r12pt);
  # ---- (F2) and the singular row at t = 0
  chk("F2", "r01(z;0) = c/(c+z), z^$_", $r01->[$_][0], R((-1) ** $_) / $c ** $_) for 0..$MZ;
  # ---- 10.6 (a), (b): v_2(z;0) = h^-2 (1, x/c), [t^1] r02 = 0
  chk("10.6(a)", "r02[m](0) = [z^m] h^-2, m=$_", $r02->[$_][0], $hm2->[$_]) for 0..$MZ;
  chk("10.6(a)", "point frame r02pt[m](0) = [z^m] h^-1, m=$_", $r02pt[$_][0], $hinv->[$_]) for 0..$MZ;
  chk("10.6(b)", "r12[m](0) = (x/c) r02[m](0), m=$_", $r12->[$_][0], $x / $c * $r02->[$_][0]) for 0..$MZ;
  chk("10.6'", "[t^1] r02[m] = 0, m=$_", $r02->[$_][1], R(0)) for 0..$MZ;
  # ---- Row 2 of the QDE (V5) at both idempotents
  for my $m (1..$MZ) {
    my $l2 = sup(ssub(sder($r12->[$m-1]), smul($dl2, $r12->[$m-1])), 1); my $rr2 = smul($tL1, ssub($r12->[$m], smul($L2, $r02->[$m])));
    chk("Row 2 at U2 (Lemma 10.5)", "m=$m t^$_", $l2->[$_], $rr2->[$_]) for 0..$NT - 2;
    my $l1 = ssub(sder($T1->[$m-1]), smul($dl1, $T1->[$m-1])); my $rr1 = smul($L2, ssub($T1->[$m], smul($tL1, $r01->[$m])));
    chk("Row 2 at U1", "m=$m t^$_", $l1->[$_], $rr1->[$_]) for 0..$NT - 3; }
  # ---- pairings. <a,b> = x^2 (a0 b1 + a1 b0) - c x a1 b1
  my $rr01 = [ map { smul($rho, $r01->[$_]) } 0..$MZ ];      # rho r01
  my $rT1  = [ map { smul($rho, $T1->[$_]) } 0..$MZ ];       # t rho r11
  my $trr01 = [ map { sup($rr01->[$_], 1) } 0..$MZ ];        # t rho r01
  # t <rho v1[p], v2[q]>
  my $P12t = sub { my ($p, $q) = @_; return ssub(ssc($x * $x, sadd(smul($trr01->[$p], $r12->[$q]), smul($rT1->[$p], $r02->[$q]))), ssc($c * $x, smul($rT1->[$p], $r12->[$q]))); };
  # t^2 rho^2 <v1[p], v1[q]>
  my $P11tt = sub { my ($p, $q) = @_; return ssub(ssc($x * $x, sadd(smul($trr01->[$p], $rT1->[$q]), smul($rT1->[$p], $trr01->[$q]))), ssc($c * $x, smul($rT1->[$p], $rT1->[$q]))); };
  my $P22 = sub { my ($p, $q) = @_; return ssub(ssc($x * $x, sadd(smul($r02->[$p], $r12->[$q]), smul($r12->[$p], $r02->[$q]))), ssc($c * $x, smul($r12->[$p], $r12->[$q]))); };
  # symplecticity / orthogonality (V6), K_1 = 1
  for my $n (0..$MZ) { my ($o, $k1, $k2) = (sz(), sz(), sz());
    for my $q (0..$n) { my $sg = R((-1) ** $q); $o = sadd($o, ssc($sg, $P12t->($n - $q, $q))); $k1 = sadd($k1, ssc($sg, $P11tt->($n - $q, $q))); $k2 = sadd($k2, ssc($sg, $P22->($n - $q, $q))); }
    chk("orthogonality (I5)", "t<rho v1, v2>, w^$n, t^$_", $o->[$_], R(0)) for 0..$NT - 1;
    chk("K_1 = 1", "w^$n t^$_", $k1->[$_], ($n == 0 && $_ == 0) ? -$x * $c ** 3 : R(0)) for 0..$NT - 1;
    chk("K_2 = 1", "w^$n t^$_", $k2->[$_], $n == 0 ? $D2->[$_] : R(0)) for 0..$NT - 2; }
  # ---- edges from EdN: E(m,n) = Sum_k (-1)^k P_(m+1+k, n-k), P_pq = -(-1)^(p+q) <v_i[p], v_j[q]>
  my (%tE12, %tE21, %E22, %ttE11);
  for my $m (0..$MZ - 1) { for my $n (0..$MZ - 1 - $m) { my ($s12, $s21, $s22, $s11) = (sz(), sz(), sz(), sz());
    for my $k (0..$n) { my ($p, $q) = ($m + 1 + $k, $n - $k); my $sg = R((-1) ** $k) * R(-((-1) ** ($p + $q)));
      $s12 = sadd($s12, ssc($sg, $P12t->($p, $q))); $s21 = sadd($s21, ssc($sg, $P12t->($q, $p)));
      $s22 = sadd($s22, ssc($sg, $P22->($p, $q))); $s11 = sadd($s11, ssc($sg, $P11tt->($p, $q))); }
    $tE12{"$m,$n"} = $s12; $tE21{"$m,$n"} = $s21; $E22{"$m,$n"} = $s22; $ttE11{"$m,$n"} = $s11; } }
  # regularity (F3): t rho E12, t rho E21 vanish at t = 0; t^2 rho^2 E11 = O(t^2)
  for my $key (sort keys %tE12) {
    chk("F3: [t^-1] rho E12 = 0 (Prop 9.3(i))", $key, $tE12{$key}->[0], R(0));
    chk("F3: [t^-1] rho E21 = 0 (needs 10.6(b))", $key, $tE21{$key}->[0], R(0));
    chk("F3: t^2 rho^2 E11 = O(t^2)", "$key t^0", $ttE11{$key}->[0], R(0)); chk("F3: t^2 rho^2 E11 = O(t^2)", "$key t^1", $ttE11{$key}->[1], R(0)); }
  my %rE12 = map { $_ => sdown($tE12{$_}, 1, "E12 $_") } keys %tE12;         # rho E12, exact to t^(NT-2)
  my %rE21 = map { $_ => sdown($tE21{$_}, 1, "E21 $_") } keys %tE21;
  my %r2E11 = map { $_ => sdown($ttE11{$_}, 2, "E11 $_") } keys %ttE11;      # rho^2 E11, exact to t^(NT-3)
  # transposition (Prop 9.5)
  for my $m (0..$MZ - 1) { for my $n (0..$MZ - 1 - $m) { chk("E21(b,a) = E12(a,b) (Prop 9.5)", "($m,$n) t^$_", $rE21{"$m,$n"}->[$_], $rE12{"$n,$m"}->[$_]) for 0..$NT - 2; } }
  # E11 symmetry
  for my $m (0..$MZ - 1) { for my $n (0..$MZ - 1 - $m) { chk("E11 symmetric", "($m,$n) t^$_", $r2E11{"$m,$n"}->[$_], $r2E11{"$n,$m"}->[$_]) for 0..$NT - 3; } }
  # ---- E22 at t = 0: S4 10.6(c) versus conifold_regularity_note.md F3
  my ($agreeS4, $agreeNote) = (0, 0); my $tot22 = 0;
  for my $m (0..$MZ - 1) { for my $n (0..$MZ - 1 - $m) { my $pred = R(0);
      for my $k (0..$n) { $pred += R((-1) ** $k) * (-$h2->[$m + 1 + $k] * $h2->[$n - $k]); } $pred *= $x ** 3 / $c;
      my $note = -$x ** 3 / 2 / $c ** ($m + $n + 2); my $v = $E22{"$m,$n"}->[0]; $tot22++;
      chk("E22(t=0) = (x^3/c)(1 - h^2 h^2)/(z+w) (10.6(c))", "($m,$n)", $v, $pred); $agreeNote++ if $v == $note; } }
  print "  E22(0,0)|t=0 = $E22{'0,0'}->[0];  S4: ", -$x * $x * (1 + $x + $x * $x) / (6 * $c * $c), ";  note F3: ", -$x ** 3 / (2 * $c * $c), "  (note F3 agrees at $agreeNote of $tot22 entries)\n";
  # ---- closed forms (S4.4) for comparison: A_j, l, S_j^(k)
  my $dphi = sint(ssub(sc(1), ssc($x, $sst)));
  my $omega = smul($tot, sexp(ssc(R(1) / $c, $dphi)));
  my @A; for my $j (1..$NT + 1) { my $C = R(1); $C *= ($j - $_ * $c) ** 2 for 1..$j - 1; $C /= fact($j) * fact($j - 1); $C /= $c ** (2 * $j - 4);
    my $om = sc(1); $om = smul($om, $omega) for 1..$j; $A[$j] = sup(ssc(R((-1) ** $j) * $C, smul($om, scompose(hyp21($j * $x / $c, $j * $x / $c, R($j)), $tau))), $j - 1); }
  my $l = sdown(ssub($tL1, sc($c)), 1, "l");                                   # L1 - c/t, exact to NT-1
  my $aj = sub { my ($j, $k) = @_; return ($k < 0 || $j > $NT + 1) ? R(0) : $A[$j]->[$k]; };
  my $Sj = sub { my ($j, $k) = @_; my $s1 = R($k + 2 - $j) * $aj->($j, $k + 1); $s1 -= R($j) / $c * $l->[$_] * $aj->($j, $k - $_) for 0..$k; return (R(0 - $j) / $c * $aj->($j, $k), $s1); };
  # (9.1) against the Wick rows:  [t^k][psi^a] sigma = (-1)^a [t^k] rho v1[a] = Sum_j S_j^(k) j^a/c^(a+1)
  for my $k (0..$NT - 3) { for my $a (0..$MZ) { my ($p0, $p1) = (R(0), R(0));
      for my $j (1..$k + 1) { my ($u0, $u1) = $Sj->($j, $k); $p0 += $u0 * R($j) ** $a / $c ** ($a + 1); $p1 += $u1 * R($j) ** $a / $c ** ($a + 1); }
      chk("(9.1) sigma_k leaf form vs Wick rows", "k=$k a=$a comp0", R((-1) ** $a) * $rr01->[$a][$k], $p0);
      chk("(9.1) sigma_k leaf form vs Wick rows", "k=$k a=$a comp1", R((-1) ** $a) * $rT1->[$a][$k + 1], $p1); } }
  for my $a (0..$MZ) { chk("(9.1) sigma_-1 = (0, c)", "a=$a", R((-1) ** $a) * $rT1->[$a][0], $a == 0 ? $c : R(0)); }
  # ---- E12 legs: direct (EdN) versus (9.8) versus the Mathematica fits
  for my $b (0..$MZ - 1) { for my $k (0..$NT - 3) {
      my @X; for my $j (1..$k + 1) { my $Xj = R(0);
        for my $k1 ($j - 1 .. $k) { my ($s0, $s1) = $Sj->($j, $k1); my $k2 = $k - $k1; my ($y0, $y1) = (R(0), R(0));
          for my $q (0..$b) { my $f = R((-1) ** $q) * R(-$j) ** ($b - $q) / $c ** ($b - $q + 1); $y0 += $f * $r02->[$q][$k2]; $y1 += $f * $r12->[$q][$k2]; }
          $Xj += $x * $x * ($s0 * $y1 + $s1 * $y0) - $c * $x * $s1 * $y1; }
        $X[$j] = $Xj; }
      for my $a (0..$MZ - 1 - $b) { my $pred = R(0); for my $j (1..$k + 1) { $pred -= R($j) * R($j) ** $a / $c ** ($a + 1) * $X[$j]; }
        chk("(9.8) E12 legs vs EdN", "b=$b k=$k a=$a", $rE12{"$a,$b"}->[$k], $pred); }
      if (exists $leg{"2,$b,$k"}) { my @e = @{$leg{"2,$b,$k"}};
        for my $a (0..$MZ - 1 - $b) { my $v = R(0); $v += evalMma($e[$_], $x) * R($_) ** $a / $c ** ($a + 1) for 0..$#e;
          chk("E12 legs: EdN vs Mathematica fit", "b=$b k=$k a=$a", $rE12{"$a,$b"}->[$k], $v); } }
      if (exists $leg{"1,$b,$k"}) { my @e = @{$leg{"1,$b,$k"}}; my $rinv = sinv($rho);
        for my $a (0..$MZ - 1 - $b) { my $v = R(0); $v += evalMma($e[$_], $x) * R($_) ** $a / $c ** ($a + 1) for 0..$#e;
          my $direct = smul($rinv, $r2E11{"$a,$b"})->[$k];                          # rho f_{1,b} = rho^-1 [w^b] rho^2 E11
          chk("E11 legs: EdN vs Mathematica fit", "b=$b k=$k a=$a", $direct, $v) if $k <= $NT - 3; } } } }
  # ---- loops: rho^2 E11 versus CF.4  (e_il = (il/(c^2 t^2)) Int Theta B_i' B_l')
  my $Th = ssc(-$x ** 3, smul($tau, smul($sst, $sst)));
  my %Bd; $Bd{$_} = sder(sup($A[$_], 1)) for 1..$NT;
  for my $k (0..$NT - 3) { for my $m (0..$MZ - 1) { for my $n (0..$MZ - 1 - $m) { my $pred = R(0);
        for my $i (1..$k + 1) { for my $ll (1..$k + 1) { next if ($i - 1) + ($ll - 1) > $k;
            my $e = ssc(R($i * $ll) / ($c * $c), sdown(sint(smul($Th, smul($Bd{$i}, $Bd{$ll}))), 2, "e"))->[$k];
            $pred += $e * R($i) ** $m * R($ll) ** $n / $c ** ($m + $n + 2); } }
        chk("CF.4 loops vs EdN", "k=$k ($m,$n)", $r2E11{"$m,$n"}->[$k], $pred); } } }
  printf "  done x = %s (%d s)\n", $x, time - $t0;
}
print "==================== summary (comparisons / mismatches)\n";
printf "  %-58s %5d  %d\n", $_, $nC{$_}, $nB{$_} // 0 for sort keys %nC;
my $bad = 0; $bad += $_ for values %nB; my $tot = 0; $tot += $_ for values %nC;
print "  total comparisons: $tot (x = @XS)\n"; print $bad ? "  TOTAL MISMATCHES: $bad\n" : "  all comparisons exact\n";
