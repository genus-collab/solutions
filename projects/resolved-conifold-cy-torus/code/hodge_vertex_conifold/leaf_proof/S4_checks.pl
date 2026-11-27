#!/usr/bin/perl
# S4_checks.pl -- exact rational spot checks (Perl + Math::BigRat, no Mathematica) for section S4 of the leaf proof
# (leaf_proof/S4_legs_closed.md).  At a rational value of x it computes, as exact truncated power series in t:
#   tau = 1 - e^-t, s* (root of x tau s^2 - c s + 1 = 0, s*(0) = 1/c), U1 = 1 - tau s*, dphi (three expressions),
#   omega = (tau/t) e^(dphi/c), the leaf coefficients A_j(t) = (-1)^j C_j t^(j-1) omega^j 2F1(jx/c, jx/c; j; tau),
#   rho (three expressions), Theta (four expressions), Delta_2 (two expressions), the loop coefficients
#   e_il(t) = (i l/(c^2 t^2)) Int_0^t Theta (s A_i)' (s A_l)' ds, the E11 and E12 legs (t^0 and top coefficients),
# and compares with (a) internal identities (sum rule Sum_j A_j = -c^2, partial fractions of [t^k](z G) evaluated directly
# from the KEY FORMULA at rational z) and (b) the Mathematica tables hvc_leaf_structure_out.txt,
# hvc_leaf_structure_out_kt5.txt, hvc_leaf_proof_check_out.txt (parsed and evaluated at the same x).
# Usage: perl S4_checks.pl [x ...]      (default: 3/7 2/5 5/3 7/2)
use strict; use warnings; use Math::BigRat; use FindBin;

my $ROOT = "$FindBin::Bin/..";
our $M = 11;                                   # series are kept to t^M
my @XS = @ARGV ? @ARGV : ("3/7", "2/5", "5/3", "7/2");

# ---------------------------------------------------------------- truncated power series (array refs of BigRat, indices 0..M)
sub R { my $v = shift; return ref($v) ? $v->copy : Math::BigRat->new("$v"); }
sub sz { return [ map { Math::BigRat->new(0) } 0..$M ]; }
sub sc { my $c = shift; my $s = sz(); $s->[0] = R($c); return $s; }
sub sadd { my ($a, $b) = @_; return [ map { $a->[$_] + $b->[$_] } 0..$M ]; }
sub ssub { my ($a, $b) = @_; return [ map { $a->[$_] - $b->[$_] } 0..$M ]; }
sub ssc  { my ($k, $a) = @_; my $kk = R($k); return [ map { $a->[$_] * $kk } 0..$M ]; }
sub smul { my ($a, $b) = @_; my $r = sz();
  for my $i (0..$M) { next if $a->[$i]->is_zero; for my $j (0..$M-$i) { next if $b->[$j]->is_zero; $r->[$i+$j] += $a->[$i] * $b->[$j]; } }
  return $r; }
sub sinv { my $a = shift; die "sinv: zero constant term" if $a->[0]->is_zero; my $r = sz(); $r->[0] = R(1) / $a->[0];
  for my $n (1..$M) { my $s = R(0); for my $k (1..$n) { $s += $a->[$k] * $r->[$n-$k]; } $r->[$n] = -$s / $a->[0]; } return $r; }
sub sdiv { my ($a, $b) = @_; return smul($a, sinv($b)); }
sub sder { my $a = shift; my $r = sz(); for my $n (1..$M) { $r->[$n-1] = $a->[$n] * $n; } return $r; }        # exact to t^(M-1)
sub sint { my $a = shift; my $r = sz(); for my $n (0..$M-1) { $r->[$n+1] = $a->[$n] / ($n+1); } return $r; }  # zero constant
sub sexp { my $a = shift; die "sexp: nonzero constant term" unless $a->[0]->is_zero; my $r = sz(); $r->[0] = R(1);
  for my $n (1..$M) { my $s = R(0); for my $k (1..$n) { $s += $a->[$k] * $r->[$n-$k] * $k; } $r->[$n] = $s / $n; } return $r; }
sub slog { my $a = shift; die "slog: constant term not 1" unless $a->[0] == 1; return sint(smul(sder($a), sinv($a))); }
sub spow { my ($a, $p) = @_; die "spow: constant term not 1" unless $a->[0] == 1; return sexp(ssc($p, slog($a))); }
sub sup  { my ($a, $m) = @_; my $r = sz(); for my $n ($m..$M) { $r->[$n] = $a->[$n-$m]->copy; } return $r; }    # times t^m
sub sdown { my ($a, $m) = @_; for my $n (0..$m-1) { die "sdown: nonzero low coefficient" unless $a->[$n]->is_zero; }
  my $r = sz(); for my $n (0..$M-$m) { $r->[$n] = $a->[$n+$m]->copy; } return $r; }                            # divide by t^m (top m lost)
sub scompose { my ($f, $g) = @_; die "compose: g(0) != 0" unless $g->[0]->is_zero; my $r = sz(); my $p = sc(1);
  for my $n (0..$M) { $r = sadd($r, ssc($f->[$n], $p)); $p = smul($p, $g); } return $r; }
sub hyp21 { my ($a, $b, $g) = @_; my @f = (R(1));
  for my $n (1..$M) { push @f, $f[-1] * ($a + $n - 1) * ($b + $n - 1) / (($g + $n - 1) * $n); } return \@f; }
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 1..$n; return $f; }

# ---------------------------------------------------------------- parsing of Mathematica InputForm (evaluated with bigrat)
sub splitTop { my $s = shift; my @out; my $d = 0; my $cur = "";
  for my $ch (split //, $s) { if ($ch eq '(' || $ch eq '{' || $ch eq '[') { $d++; } elsif ($ch eq ')' || $ch eq '}' || $ch eq ']') { $d--; }
    if ($ch eq ',' && $d == 0) { push @out, $cur; $cur = ""; } else { $cur .= $ch; } }
  push @out, $cur if $cur =~ /\S/; s/^\s+|\s+$//g for @out; return @out; }
sub evalMma { my ($expr, $xv, $tv) = @_; my $p = $expr; $p =~ s/\^/**/g; $p =~ s/\bx\b/\$x/g; $p =~ s/\bt\b/\$t/g;
  my $x = $xv; my $t = defined $tv ? $tv : 0; my $v; { use bigrat; $v = eval $p; } die "eval failed ($@): $p" if $@ || !defined $v;
  return Math::BigRat->new("$v"); }
sub solveLin { my ($A, $b) = @_; my $n = scalar @$b; my @a = map { [ (map { R($_) } @{$A->[$_]}), R($b->[$_]) ] } 0..$n-1;
  for my $i (0..$n-1) { my $p = $i; $p++ while $a[$p][$i]->is_zero; @a[$i, $p] = @a[$p, $i];
    for my $r (0..$n-1) { next if $r == $i; my $f = $a[$r][$i] / $a[$i][$i]; next if $f->is_zero; $a[$r][$_] -= $f * $a[$i][$_] for $i..$n; } }
  return [ map { $a[$_][$n] / $a[$_][$_] } 0..$n-1 ]; }
sub tPolyCoeffs { my ($expr, $xv, $deg) = @_; my (@A, @b);           # coefficients of a polynomial in t (degree <= deg) by interpolation
  for my $i (1..$deg+1) { my $tv = R($i); push @A, [ map { $tv ** $_ } 0..$deg ]; push @b, evalMma($expr, $xv, $tv); }
  return solveLin(\@A, \@b); }
sub slurp { my $f = shift; open(my $fh, "<", $f) or die "cannot open $f"; local $/; my $s = <$fh>; close $fh; return $s; }

my $outA = slurp("$ROOT/hvc_leaf_structure_out.txt");
my $outB = slurp("$ROOT/hvc_leaf_structure_out_kt5.txt");
my $outC = slurp("$ROOT/hvc_leaf_proof_check_out.txt");
my (%dil, %leg, %loop, @rhoTab, $dphiTab, $thetaTab);
for my $txt ($outA, $outB) { for my $line (split /\n/, $txt) {
  if ($line =~ /^\(i\) rho q_t - q_0, t\^(\d+): coefficients of lam_0\.\.lam_\d+ = \{(.*)\}; checks/) { $dil{$1} = [ splitTop($2) ]; }
  elsif ($line =~ /^\(ii\) leg j = (\d+), b = (\d+), t\^(\d+): lam_0\.\.lam_\d+ = \{(.*)\}; checks/) { $leg{"$1,$2,$3"} = [ splitTop($4) ]; }
  elsif ($line =~ /^\(iii'\) rho\^2 E11, t\^(\d+), restricted ansatz \(no lam_0, D <= k\): \{(.*)\}; extra equations/) {
    my $k = $1; for my $p (splitTop($2)) { $p =~ /^\{(\d+), (\d+)\} -> (.*)$/ or die "loop entry: $p"; $loop{"$k,$1,$2"} = $3; } }
  elsif ($line =~ /^rho = \{(.*)\}/) { my @r = splitTop($1); @rhoTab = @r if @r > @rhoTab; } } }
($dphiTab) = $outC =~ /^dphi\(t\) = (.*)$/m; ($thetaTab) = $outC =~ /^Theta\(t\) = (.*)$/m;

my @B2 = (undef, R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"));   # B_2k, k = 1..6

my (%nCmp, %nBad);
sub cmpv { my ($cat, $label, $a, $b) = @_; $nCmp{$cat}++; if ($a != $b) { $nBad{$cat}++; print "  MISMATCH [$cat] $label: $a vs $b\n"; } }

for my $xs (@XS) {
  my $x = R($xs); my $c = 1 + $x; print "==================== x = $x (c = $c)\n";
  # tau and tau/t
  my $tau = sz(); my $tot = sz();
  for my $n (0..$M) { $tot->[$n] = R((-1) ** $n) / fact($n + 1); $tau->[$n + 1] = $tot->[$n]->copy if $n + 1 <= $M; }
  my $tto = sinv($tot);                                                   # t/tau
  # s*: s = (1 + x tau s^2)/c
  my $s = sc(R(1) / $c); for (0..$M + 1) { $s = ssc(R(1) / $c, sadd(sc(1), ssc($x, smul($tau, smul($s, $s))))); }
  my $crit = sadd(ssub(ssc($x, smul($tau, smul($s, $s))), ssc($c, $s)), sc(1)); cmpv("algebra", "critical equation t^$_", $crit->[$_], R(0)) for 0..$M;
  my $ts = smul($tau, $s); my $U1 = ssub(sc(1), $ts);
  # dphi: (a) Int (1 - x s*) dt, (b) c log(c s*) + (x - 1) log(1 - tau s*), (c) definition phi_tau(s*) - phi_0(1/c)
  my $dphiA = sint(ssub(sc(1), ssc($x, $s)));
  my $dphiB = sadd(ssc($c, slog(ssc($c, $s))), ssc($x - 1, slog($U1)));
  my $dphiC = ssub(sadd(slog(ssc($c, $s)), ssc($x, slog(ssc($c / $x, ssub(sc(1), $s))))), slog($U1));
  for my $k (0..$M) { cmpv("dphi", "A vs B t^$k", $dphiA->[$k], $dphiB->[$k]); cmpv("dphi", "A vs C t^$k", $dphiA->[$k], $dphiC->[$k]); }
  my $dT = tPolyCoeffs($dphiTab, $x, 7); cmpv("dphi", "table t^$_", $dphiA->[$_], $dT->[$_]) for 0..7;
  # omega and the leaf coefficients A_j
  my $omega = smul($tot, sexp(ssc(R(1) / $c, $dphiA)));
  my (@A, @Cj);
  for my $j (1..$M + 1) {
    my $C = R(1); $C *= ($j - $_ * $c) ** 2 for 1..$j - 1; $C /= fact($j) * fact($j - 1); $C /= $c ** (2 * $j - 4); $Cj[$j] = $C;
    my $om = sc(1); $om = smul($om, $omega) for 1..$j;
    my $F = scompose(hyp21($j * $x / $c, $j * $x / $c, R($j)), $tau);
    $A[$j] = sup(ssc(R((-1) ** $j) * $C, smul($om, $F)), $j - 1); }
  # (1) sum rule Sum_j A_j = -c^2
  for my $k (0..$M) { my $sum = R(0); $sum += $A[$_][$k] for 1..$M + 1; cmpv("sum rule", "t^$k", $sum, $k == 0 ? -$c * $c : R(0)); }
  # (2) partial fractions: [t^k](z G) from the KEY FORMULA at rational z  =  c delta_k0 + Sum_j A_j[k]/(c + j z)
  for my $zs ("1/3", "2", "-1/5", "7/2", "5", "3/11", "-9/4", "1/17", "13/3", "-2/7", "11", "4/9", "-6", "8/5") {
    my $z = R($zs); next if grep { ($c + $_ * $z)->is_zero } 1..$M + 1;    # z must avoid the poles -c/j
    my @Fn = (R(1));
    for my $n (1..$M) { my $v = R(1); $v *= (1 + $_ * $z) ** 2 for 1..$n; my $d = fact($n) * $z ** $n; $d *= ($c + $_ * $z) for 2..$n + 1; push @Fn, $v / $d; }
    my $zG = ssc($c * $z / ($c + $z), smul($tot, smul(scompose(\@Fn, $tau), sexp(ssc(R(-1) / $z, $dphiA)))));    # (not -R(1): that is the file test -R)
    for my $k (0..$M) { my $rhs = $k == 0 ? $c->copy : R(0); $rhs += $A[$_][$k] / ($c + $_ * $z) for 1..$M + 1; cmpv("partial fractions", "z=$z t^$k", $zG->[$k], $rhs); } }
  # (3) dilaton tables: delta_k = Sum_j A_j[k] lam_j (lam_0 coefficient 0)
  for my $k (sort { $a <=> $b } keys %dil) { my @e = @{$dil{$k}};
    cmpv("delta_k table", "k=$k lam_0", evalMma($e[0], $x), R(0));
    for my $j (1..$#e) { cmpv("delta_k table", "k=$k lam_$j", evalMma($e[$j], $x), $A[$j][$k] - (($k == 0 && $j == 1) ? -$c * $c : 0)); } }
  # (4) top and subleading coefficients, all k <= M
  for my $k (0..$M) { cmpv("top coefficient", "a_{k+1}^(k), k=$k", $A[$k + 1][$k], R((-1) ** ($k + 1)) * $Cj[$k + 1]); }
  for my $k (1..$M) { cmpv("subleading coefficient", "a_k^(k), k=$k", $A[$k][$k], R((-1) ** $k) * $k * $Cj[$k] * (1 - $x) ** 2 / (2 * $c * $c)); }
  # (5) rho: closed form, from Delta_1 (V2 formula), from alg; table
  my $u = ssc(4 * $x / ($c * $c), $tau);
  my $rhoClosed = smul($tot, ssc(2, sinv(smul(spow(ssub(sc(1), $u), R("1/4")), sadd(sc(1), spow(ssub(sc(1), $u), R("1/2")))))));
  my $num = ssub(ssc($c, $tto), ssc(2 * $x, sup($s, 1)));                 # t^2 (L1 - L2)/t
  my $den = sdown(ssc($x, ssub(sc(1), ssc($c, $s))), 1);                 # (x - c L2)/t, exact to t^(M-1)
  my $d1t2 = ssc($x ** 3, sdiv($num, $den));                             # Delta_1 t^2, exact to t^(M-1)
  my $rhoV2 = spow(ssc(R(1) / (-$x * $c ** 3), $d1t2), R("-1/2"));
  my $phi2 = sadd(ssub(ssc(-1, sinv(smul($s, $s))), ssc($x, sinv(smul(ssub(sc(1), $s), ssub(sc(1), $s))))), smul(smul($tau, $tau), sinv(smul($U1, $U1))));
  my $alg = smul($U1, spow(ssc(-$x / $c ** 3, $phi2), R("1/2")));
  my $rhoAlg = sdiv($tot, $alg);
  for my $k (0..$M - 1) { cmpv("rho", "closed vs V2 t^$k", $rhoClosed->[$k], $rhoV2->[$k]); cmpv("rho", "closed vs alg t^$k", $rhoClosed->[$k], $rhoAlg->[$k]); }
  cmpv("rho", "table t^$_", $rhoClosed->[$_], evalMma($rhoTab[$_], $x)) for 0..$#rhoTab;
  # (6) Theta: -x^3 tau s*^2, x c dphi' - x, -4 x^3 tau/(c + delta)^2, V8 definition; table
  my $thA = ssc(-$x ** 3, smul($tau, smul($s, $s)));
  my $thB = ssub(ssc($x * $c, ssub(sc(1), ssc($x, $s))), sc($x));
  my $delta = ssc($c, spow(ssub(sc(1), $u), R("1/2")));
  my $thC = ssc(-4 * $x ** 3, sdiv($tau, smul(sadd(sc($c), $delta), sadd(sc($c), $delta))));
  my $q = ssub(ssc($c, $s), sc(1)); my $q2t = smul(sdown(smul($q, $q), 1), $tto);     # (c s* - 1)^2/tau, exact to t^(M-1)
  my $thD = ssc($x, smul($q2t, sadd(ssub(ssc($x, $tau), sc($c * $c)), ssc($c * $x, $ts))));
  for my $k (0..$M - 1) { cmpv("Theta", "A vs B t^$k", $thA->[$k], $thB->[$k]); cmpv("Theta", "A vs C t^$k", $thA->[$k], $thC->[$k]); cmpv("Theta", "A vs V8 t^$k", $thA->[$k], $thD->[$k]); }
  my $thT = tPolyCoeffs($thetaTab, $x, 7); cmpv("Theta", "table t^$_", $thA->[$_], $thT->[$_]) for 0..7;
  # (7) Delta_2 = x^3 (L2 - L1)/(x - c L1) = x^3 delta s*^2; Delta_2(0) = x^3/c
  my $d2a = ssc($x ** 3, sdiv(ssub(ssc(2 * $x, $ts), sc($c)), sadd(ssub(ssc($x, $tau), sc($c * $c)), ssc($c * $x, $ts))));
  my $d2b = ssc($x ** 3, smul($delta, smul($s, $s)));
  cmpv("Delta_2", "V2 formula vs x^3 delta s*^2 t^$_", $d2a->[$_], $d2b->[$_]) for 0..$M; cmpv("Delta_2", "value at 0", $d2a->[0], $x ** 3 / $c);
  # (8) second component of sigma: S_j = (-(j/c) A_j, A_j' - ((j-1)/t) A_j - (j l/c) A_j); l = L1 - c/t; S_j = O(t^(j-1)); sigma_0 = (c, x) lam_1
  my $l = ssub(ssc($c, sdown(ssub($tto, sc(1)), 1)), ssc($x, $s));     # c (1/tau - 1/t) - x s*, exact to t^(M-1)
  cmpv("sigma", "l(0) = (1 + x^2)/(2c)", $l->[0], (1 + $x * $x) / (2 * $c));
  for my $j (1..6) {   # A_j/t is a power series for j >= 2 (A_j = O(t^(j-1))); for j = 1 the term (j-1) A_j/t is absent
    my $S1 = $j >= 2 ? ssub(sder($A[$j]), ssc($j - 1, sdown($A[$j], 1))) : sder($A[$j]);
    $S1 = ssub($S1, ssc(R($j) / $c, smul($l, $A[$j])));
    cmpv("sigma", "S_$j second component vanishes below t^" . ($j - 1) . " (t^$_)", $S1->[$_], R(0)) for 0..$j - 2;
    if ($j == 1) { cmpv("sigma", "S_1(0) = (c, x)", -$A[1][0] / $c, $c); cmpv("sigma", "S_1(0) second", $S1->[0], $x); } }
  # (9) loops: e_il(t) = (i l/(c^2 t^2)) Int_0^t Theta B_i' B_l' (B_j = t A_j); table kt5; corner and anti-diagonal
  my (%e, %Bd); $Bd{$_} = sder(sup($A[$_], 1)) for 1..8;
  for my $i (1..7) { for my $ll ($i..7) { $e{"$i,$ll"} = ssc(R($i * $ll) / ($c * $c), sdown(sint(smul($thA, smul($Bd{$i}, $Bd{$ll}))), 2)); $e{"$ll,$i"} = $e{"$i,$ll"}; } }
  for my $key (sort keys %loop) { my ($k, $i, $ll) = split /,/, $key; cmpv("E11 loops table", "t^$k ($i,$ll)", $e{"$i,$ll"}[$k], evalMma($loop{$key}, $x)); }
  for my $k (0..5) { for my $i (1..$k + 1) { for my $ll ($i..$k + 1) { next if ($i - 1) + ($ll - 1) <= $k; cmpv("E11 support", "t^$k ($i,$ll) outside support", $e{"$i,$ll"}[$k], R(0)); } } }
  for my $k (0..6) { cmpv("E11 corner", "k=$k", $e{"1," . ($k + 1)}[$k], R(($k + 1) ** 2) / ($k + 2) * $x ** 3 / ($c * $c) * $A[$k + 1][$k]); }
  for my $i (1..5) { for my $ll ($i..5) { my $k = $i + $ll - 2; next if $k > $M - 3;
    cmpv("E11 anti-diagonal", "($i,$ll) at t^$k", $e{"$i,$ll"}[$k], -$x ** 3 * $i * $i * $ll * $ll * $A[$i][$i - 1] * $A[$ll][$ll - 1] / ($c ** 4 * ($i + $ll))); } }
  # (10) E11 legs: lam_i coefficient of [t^k] rho f_{1,b} = Sum_{k1+k2=k} (1/rho)_{k1} Sum_l e_il[k2] l^b/c^(b+1)
  my $rinv = sinv($rhoClosed);
  for my $key (sort keys %leg) { my ($jj, $b, $k) = split /,/, $key; next unless $jj == 1; my @ent = @{$leg{$key}};
    cmpv("E11 legs table", "b=$b t^$k lam_0", evalMma($ent[0], $x), R(0));
    for my $i (1..$#ent) { my $v = R(0);
      for my $k2 (0..$k) { for my $ll (1..7) { $v += $rinv->[$k - $k2] * $e{"$i,$ll"}[$k2] * R($ll) ** $b / $c ** ($b + 1); } }
      cmpv("E11 legs table", "b=$b t^$k lam_$i", evalMma($ent[$i], $x), $v); } }
  # (11) E12 legs: t^0 value chi_b = -x^3 [w^b] h_w(w)^2/(c + w); top coefficient (k+1)^2 (x^3/c^2) a_{k+1}^(k) [w^b] h_w(w)^2/(c + (k+1) w)
  my $ex = sz(); for my $kk (1..6) { my $m = 2 * $kk - 1; last if $m > $M; my $p = R(1) / $c ** $m - R(1) / $x ** $m - 1; $ex->[$m] = -2 * $B2[$kk] / (2 * $kk * (2 * $kk - 1)) * $p; }
  my $h2 = sexp($ex);                                                    # h_w(w)^2 as a series in w
  for my $key (sort keys %leg) { my ($jj, $b, $k) = split /,/, $key; next unless $jj == 2; my @ent = @{$leg{$key}};
    cmpv("E12 legs table", "b=$b t^$k lam_0", evalMma($ent[0], $x), R(0));
    my $jt = $k + 1; my $g = sdiv($h2, sadd(sc($c), ssc($jt, sup(sc(1), 1))));     # h_w^2/(c + (k+1) w)
    cmpv("E12 legs top", "b=$b t^$k lam_$jt", evalMma($ent[$jt], $x), R($jt * $jt) * $x ** 3 / ($c * $c) * $A[$jt][$k] * $g->[$b]);
    if ($k == 0) { cmpv("E12 legs t^0 (chi_b)", "b=$b", evalMma($ent[1], $x), -$x ** 3 * $g->[$b]); } }
}
print "==================== summary (comparisons / mismatches)\n";
printf "  %-24s %5d  %d\n", $_, $nCmp{$_}, $nBad{$_} // 0 for sort keys %nCmp;
my $bad = 0; $bad += $_ for values %nBad; my $tot = 0; $tot += $_ for values %nCmp;
print "  total comparisons: $tot (x = @XS)\n"; print $bad ? "  TOTAL MISMATCHES: $bad\n" : "  all comparisons exact\n";
