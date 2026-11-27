#!/usr/bin/perl
# Continuity lemma: exact (Math::BigRat) consistency checks. Evidence only, not proof steps.
# (A) generic alpha: Hessian of W in log-coordinates at a critical point has det = W(Y)/3,
#     and W(Y) = -Y * Lhat'(Y) on Lhat(Y)=0.   (Lhat(Y) = prod(Y+a_j) + 27 y Y^3)
# (B) T_F (alpha_j = E/3), y = c^3: the formal Laplace expansion of e^{W/z} dOmega at the
#     critical point (log coords), S = W^{-1/2}(1 + z*omega1), solves the conjugated PF
#     operator eq:pf at orders z^0, z^1, z^2.  Controls: omega1 + const (allowed, kappa(z)
#     freedom) must pass; omega1 + c (non-constant) must fail at z^2; S0 = 1/q^2 must fail at z^1.
use strict; use warnings;
use Math::BigRat;
sub R { Math::BigRat->new(shift) }

# ---------------- (A) generic alpha ----------------
print "(A) generic-alpha identities\n";
my @tests = ([R('2/3'),R('5/7'),R('-11/5'),R('3/4')], [R('1'),R('2'),R('7'),R('-1/3')], [R('7/3'),R('7/3'),R('7/3'),R('-5/2')]);
for my $tc (@tests) {
  my ($a0,$a1,$a2,$Y) = @$tc;
  my @a = ($a0,$a1,$a2);
  my $prod = R(1); $prod *= ($Y + $_) for @a;
  my $y = -$prod / (27 * $Y**3);                       # makes Y a root of Lhat
  my $e1 = $a0+$a1+$a2; my $e2 = $a0*$a1+$a0*$a2+$a1*$a2; my $e3 = $a0*$a1*$a2;
  my $W = $e1*$Y**2 + 2*$e2*$Y + 3*$e3;
  # Lhat'(Y)
  my $Lp = ($Y+$a1)*($Y+$a2) + ($Y+$a0)*($Y+$a2) + ($Y+$a0)*($Y+$a1) + 81*$y*$Y**2;
  # Hessian in x_j = log q_j : q_j delta + (P/9) J, q_j = a_j + Y, P = -3Y
  my @q = map { $_ + $Y } @a; my $P = -3*$Y;
  my @H; for my $i (0..2){ for my $j (0..2){ $H[$i][$j] = ($i==$j ? $q[$i] : R(0)) + $P/9; } }
  my $det = $H[0][0]*($H[1][1]*$H[2][2]-$H[1][2]*$H[2][1]) - $H[0][1]*($H[1][0]*$H[2][2]-$H[1][2]*$H[2][0]) + $H[0][2]*($H[1][0]*$H[2][1]-$H[1][1]*$H[2][0]);
  # criticality check: q_j - a_j + P/3 = 0 and P^3 = prod q / y
  my $crit = ($q[0]-$a0+$P/3) == 0 && ($P**3 == $q[0]*$q[1]*$q[2]/$y);
  printf "  alpha=(%s,%s,%s) Y=%s: crit=%s  det-W/3=%s  W+Y*Lhat'=%s\n", $a0,$a1,$a2,$Y, ($crit?'ok':'FAIL'), $det - $W/3, $W + $Y*$Lp;
}

# ---------------- (B) T_F, truncated Taylor series in t = c - c0 ----------------
my $N = 9;              # series length (orders 0..N-1)
my $E = R(7);           # e1 = 2, e2 = 5 (only E enters)
my $c0 = R('1/5');
sub sc { my $v = shift; my @s = (R(0)) x $N; $s[0] = R($v); return \@s; }
sub sadd { my ($a,$b) = @_; return [ map { $a->[$_] + $b->[$_] } 0..$N-1 ]; }
sub ssub { my ($a,$b) = @_; return [ map { $a->[$_] - $b->[$_] } 0..$N-1 ]; }
sub sscal { my ($k,$a) = @_; return [ map { $k * $a->[$_] } 0..$N-1 ]; }
sub smul { my ($a,$b) = @_; my @r = (R(0)) x $N;
  for my $i (0..$N-1) { next if $a->[$i] == 0; for my $j (0..$N-1-$i) { $r[$i+$j] += $a->[$i]*$b->[$j]; } } return \@r; }
sub sinv { my $a = shift; my @r = (R(0)) x $N; $r[0] = 1/$a->[0];
  for my $n (1..$N-1) { my $s = R(0); $s += $a->[$_]*$r[$n-$_] for 1..$n; $r[$n] = -$s/$a->[0]; } return \@r; }
my $cser = sc($c0); $cser->[1] = R(1);
sub sder { my $a = shift; my @r = (R(0)) x $N; $r[$_-1] = $_ * $a->[$_] for 1..$N-1; $r[$N-1] = R(0); return \@r; }
# theta = y d/dy = (c/3) d/dc ; derivative loses the top order (marked by truncation bookkeeping below)
sub theta { my $a = shift; return sscal(R('1/3'), smul($cser, sder($a))); }

my $inv1 = sinv(sadd(sc(1), sscal(R(3),$cser)));     # 1/(1+3c)
my $q  = sscal($E, smul($cser,$inv1));                # common value q_j
my $Y  = sscal(-$E/3, $inv1);                         # Lhat root, Y = q - E/3
my $P  = sscal($E, $inv1);                            # P = -3Y
my $y  = smul($cser, smul($cser,$cser));

# derivative tensors of W at the critical point (log coords): order m: q*delta + P/3^m
sub Wt { my ($m, @idx) = @_; my $s = sscal(R(1)/ (R(3)**$m), $P);
  my $alleq = 1; for (@idx) { $alleq = 0 if $_ != $idx[0]; }
  return $alleq ? sadd($s, $q) : $s; }
# K = Hess^{-1} = (1/q)(I - lam J), lam = (P/9)/(q + P/3)
my $iq = sinv($q);
my $lam = smul(sscal(R('1/9'),$P), sinv(sadd($q, sscal(R('1/3'),$P))));
my @K; for my $i (0..2){ for my $j (0..2){ $K[$i][$j] = smul($iq, ssub(sc($i==$j?1:0), $lam)); } }
# check K * Hess = I
my $kerr = 0;
for my $i (0..2){ for my $j (0..2){ my $s = sc(0); for my $k (0..2){ $s = sadd($s, smul($K[$i][$k], Wt(2,$k,$j))); }
  my $tgt = sc($i==$j?1:0); for (0..$N-1){ $kerr++ if $s->[$_] != $tgt->[$_]; } } }
print "(B) T_F at c0=$c0, E=$E: K*Hess = I exact: ", ($kerr==0?'yes':"NO ($kerr)"), "\n";

# omega1 = (1/8) W4 KK - (1/8) W3 W3 KKK(loop) - (1/12) W3 W3 KKK(cross)
my $t4 = sc(0); my $tl = sc(0); my $tx = sc(0);
for my $a (0..2){ for my $b (0..2){ for my $c (0..2){ for my $d (0..2){
  $t4 = sadd($t4, smul(Wt(4,$a,$b,$c,$d), smul($K[$a][$b],$K[$c][$d]))); }}}}
for my $a (0..2){ for my $b (0..2){ for my $c (0..2){ my $w3a = Wt(3,$a,$b,$c);
  for my $d (0..2){ for my $e (0..2){ for my $f (0..2){ my $w3b = Wt(3,$d,$e,$f); my $ww = smul($w3a,$w3b);
    $tl = sadd($tl, smul($ww, smul($K[$a][$b], smul($K[$c][$d], $K[$e][$f]))));
    $tx = sadd($tx, smul($ww, smul($K[$a][$d], smul($K[$b][$e], $K[$c][$f])))); }}}}}}
my $om1 = sadd(sscal(R('1/8'),$t4), sadd(sscal(R('-1/8'),$tl), sscal(R('-1/12'),$tx)));
print "  omega1 Taylor coeffs in t=c-c0 (orders 0..4): ", join(", ", map {"$om1->[$_]"} 0..4), "\n";
{ my $pred = sadd(sc(R(-3)/(4*$E)), sscal(R(1)/(9*$E), sinv($cser)));   # fitted: -3/(4E) + 1/(9 E c)
  my $d = ssub($om1,$pred); print "  omega1 = -3/(4E) + 1/(9Ec) to order ", $N-1, ": ", ((grep { $_ != 0 } @$d) ? 'no' : 'yes'), "\n"; }

# z-polynomials with series coefficients: arrayref of series, index = power of z (keep 0..3)
my $ZM = 3;
sub zop { # apply (z*theta + V + k z) to f ; V series, k rational
  my ($f, $V, $k) = @_; my @g = map { sc(0) } 0..$ZM;
  for my $n (0..$ZM) { $g[$n] = sadd($g[$n], smul($V, $f->[$n]));
    if ($n+1 <= $ZM) { $g[$n+1] = sadd($g[$n+1], theta($f->[$n])); $g[$n+1] = sadd($g[$n+1], sscal($k, $f->[$n])); } }
  return \@g; }
sub zadd { my ($f,$g) = @_; return [ map { sadd($f->[$_],$g->[$_]) } 0..$ZM ]; }
sub zscal { my ($s,$f) = @_; return [ map { smul($s,$f->[$_]) } 0..$ZM ]; }
sub PFop { my $S = shift; my $a = sadd($Y, sc($E/3));
  my $f = zop(zop(zop($S,$a,R(0)),$a,R(0)),$a,R(0));
  my $g = zop(zop(zop($S,$Y,R('2/3')),$Y,R('1/3')),$Y,R(0));
  return zadd($f, zscal(sscal(R(27),$y), $g)); }
sub report { my ($name, $S) = @_; my $r = PFop($S);
  my @out; for my $n (0..2) { my $ok = 1; for my $o (0..$N-1-4) { $ok = 0 if $r->[$n][$o] != 0; } push @out, "z^$n:".($ok?'0':'NONZERO'); }
  print "  $name: ", join('  ', @out), "\n"; }
my $S0 = $iq;                       # W^{-1/2} = 1/(sqrt(E) q), constant factor dropped
my @Sz = ($S0, smul($S0,$om1), sc(0), sc(0));
report("Asym S0(1+z omega1)", \@Sz);
my @Sc1 = ($S0, smul($S0, sadd($om1, sc(R('3/11')))), sc(0), sc(0));
report("control +const (allowed)", \@Sc1);
my @Sc2 = ($S0, smul($S0, sadd($om1, $cser)), sc(0), sc(0));
report("control +c (should fail z^2)", \@Sc2);
my @Sc3 = (smul($iq,$iq), sc(0), sc(0), sc(0));
report("control S0=1/q^2 (should fail z^1)", \@Sc3);
