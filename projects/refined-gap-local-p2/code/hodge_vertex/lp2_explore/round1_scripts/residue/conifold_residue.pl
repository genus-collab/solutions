# conifold_residue.pl: header (parsing + setup) copied from leaf_proof/S3_checks.pl lines 1-124; residue-lemma control appended.
#!/usr/bin/perl
# S3_checks.pl -- exact rational spot checks (Math::BigRat via bigrat) of the formulas of section S3
# (leaf_proof/S3_leaves_loops.md) against the stored Mathematica tables:
#   hvc_leaf_structure_out_kt5.txt (rho, dilaton leaves, legs, restricted loops; t^0..t^5),
#   hvc_leaf_structure_out.txt     (full loop matrices t^0..t^3),
#   hvc_leaf_proof_check_out.txt   (Theta(t), dphi(t) to t^7).
# No Mathematica is called; x is set to the rational values 2/7, 3/5, 7/3, -5/9, the leaf variables to rational
# values away from all poles; t-series are truncated at t^10 and compared only at orders where they are complete.
# Usage (from leaf_proof/): perl S3_checks.pl .. > S3_checks_out.txt      (about 10 s)
use strict; use warnings; use bigrat;

# default: frozen copies (2026-11-27) of the three resolved-conifold tables, see ref_conifold/README.md
use FindBin ();
my $dir = shift // "$FindBin::Bin/ref_conifold";
my $N = 10;                                   # series truncation t^0..t^N (comparisons only at low order)

sub R { return Math::BigRat->new("$_[0]"); }
sub szero { return [ map { R(0) } 0..$N ]; }
sub sone { my $s = szero(); $s->[0] = R(1); return $s; }
sub sadd { my ($a,$b) = @_; return [ map { $a->[$_] + $b->[$_] } 0..$N ]; }
sub ssub { my ($a,$b) = @_; return [ map { $a->[$_] - $b->[$_] } 0..$N ]; }
sub sscal { my ($k,$a) = @_; return [ map { $k * $a->[$_] } 0..$N ]; }
sub smul { my ($a,$b) = @_; my $r = szero();
  for my $i (0..$N) { next if $a->[$i] == 0; for my $j (0..$N-$i) { $r->[$i+$j] += $a->[$i]*$b->[$j]; } } return $r; }
sub sinv { my $a = shift; die "sinv: zero constant term\n" if $a->[0] == 0; my $r = szero(); $r->[0] = 1/$a->[0];
  for my $n (1..$N) { my $s = R(0); for my $i (1..$n) { $s += $a->[$i]*$r->[$n-$i]; } $r->[$n] = -$s/$a->[0]; } return $r; }
sub sexp { my $a = shift; die "sexp: nonzero constant term\n" unless $a->[0] == 0; my $r = sone();
  for my $n (1..$N) { my $s = R(0); for my $m (1..$n) { $s += $m*$a->[$m]*$r->[$n-$m]; } $r->[$n] = $s/$n; } return $r; }
sub ssqrt1 { my $a = shift; die "ssqrt1: constant term not 1\n" unless $a->[0] == 1; my $r = sone();
  for my $n (1..$N) { my $s = $a->[$n]; for my $i (1..$n-1) { $s -= $r->[$i]*$r->[$n-$i]; } $r->[$n] = $s/2; } return $r; }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; return $r; }

# ---------- parsing of Mathematica InputForm output ----------
sub mexpr { my $s = shift; $s =~ s/\^/**/g; $s =~ s/\bx\b/\$x/g; $s =~ s/\bt\b/\$t/g; return $s; }
sub meval { my ($s, $xv, $tv) = @_; my $x = $xv; my $t = defined $tv ? $tv : R(0); my $v = eval mexpr($s); die "eval failed ($@) in: $s\n" if $@; return R($v); }
sub topsplit { my $s = shift; $s =~ s/^\s*\{//; $s =~ s/\}\s*$//; my @out; my $d = 0; my $cur = '';
  for my $ch (split //, $s) { if ($ch eq '(' || $ch eq '{') { $d++; } if ($ch eq ')' || $ch eq '}') { $d--; }
    if ($ch eq ',' && $d == 0) { push @out, $cur; $cur = ''; } else { $cur .= $ch; } } push @out, $cur if $cur =~ /\S/;
  s/^\s+|\s+$//g for @out; return @out; }
sub firstlist { my $s = shift; my $i = index($s, '{'); my $d = 0; for my $k ($i..length($s)-1) { my $ch = substr($s,$k,1);
    $d++ if $ch eq '{'; $d-- if $ch eq '}'; return substr($s, $i, $k-$i+1) if $d == 0; } die "unbalanced\n"; }
sub slurp { my $f = shift; open(my $fh, '<', $f) or die "cannot open $f\n"; my @l = <$fh>; close $fh; chomp @l; return @l; }

my @kt5 = slurp("$dir/hvc_leaf_structure_out_kt5.txt");
my @kt3 = slurp("$dir/hvc_leaf_structure_out.txt");
my @chk = slurp("$dir/hvc_leaf_proof_check_out.txt");
my ($rhoL) = grep { /^rho = / } @kt5;
my (%dilL, %legL, %loopRL, %loopFL);
for (@kt5) {
  if (/^\(i\) rho q_t - q_0, t\^(\d+): coefficients of lam_0..lam_\d+ = (.*)$/) { $dilL{$1} = firstlist($2); }
  if (/^\(ii\) leg j = (\d), b = (\d+), t\^(\d+): lam_0..lam_\d+ = (\{.*)$/) { $legL{"$1,$2,$3"} = firstlist($4); }
  if (/^\(iii'\) rho\^2 E11, t\^(\d+), restricted ansatz \(no lam_0, D <= k\): (.*)$/) { $loopRL{$1} = firstlist($2); }
}
for (@kt3) { if (/^\(iii\) rho\^2 E11, t\^(\d+): matrix .*? = (\{\{.*)$/) { $loopFL{$1} = firstlist($2); } }
my ($thetaL) = grep { /^Theta\(t\) = / } @chk; $thetaL =~ s/^Theta\(t\) = //;
my ($dphiL) = grep { /^dphi\(t\) = / } @chk; $dphiL =~ s/^dphi\(t\) = //;
my $KD = 5;                                    # data available for t^0..t^5

my ($nok, $nbad) = (0, 0); my %cat;
sub check { my ($name, $ok) = @_; my $k = $name; $k =~ s/ at (z|t|psi|\().*$//; $k =~ s/\d+/#/g; $cat{$k}++;
  if ($ok) { $nok++; } else { $nbad++; print "FAIL: $name\n"; } }
END { print "number of exact checks per type (all x values together):\n"; print "  $cat{$_}  $_\n" for sort keys %cat; }

for my $xv (R("2/7"), R("3/5"), R("7/3"), R("-5/9")) {
  my $x = $xv; my $c = 1 + $x;
  print "---- x = $x ----\n";
  # tau, u = tau/t, s*, Theta, dphi, l(t) = L_1 - c/t, alg, rho
  my $tau = szero(); $tau->[$_] = R((-1)**($_+1))/fact($_) for 1..$N;
  my $u = szero(); $u->[$_] = R((-1)**$_)/fact($_+1) for 0..$N;
  my $s = szero(); $s->[0] = 1/$c;
  for (0..$N+1) { $s = sscal(1/$c, sadd(sone(), sscal($x, smul($tau, smul($s, $s))))); }
  my $crit = sadd(ssub(sscal($x, smul($tau, smul($s,$s))), sscal($c, $s)), sone());
  check("critical equation x tau s^2 - c s + 1 = 0", (grep { $crit->[$_] != 0 } 0..$N) == 0);
  my $Theta = sscal(-$x**3, smul($tau, smul($s, $s)));
  my $dphidot = ssub(sone(), sscal($x, $s));
  my $dphi = szero(); $dphi->[$_] = $dphidot->[$_-1]/$_ for 1..$N;
  # d(dphi)/dtau = s/(1 - tau s) versus d(dphi)/dt = 1 - x s: (1 - tau) s/(1 - tau s) = 1 - x s
  my $lhs = smul(ssub(sone(), $tau), smul($s, sinv(ssub(sone(), smul($tau, $s)))));
  check("(1-tau) s/(1 - tau s) = 1 - x s", (grep { $lhs->[$_] != $dphidot->[$_] } 0..$N) == 0);
  # Theta = x c d(dphi)/dt - x
  my $th2 = ssub(sscal($x*$c, $dphidot), sscal($x, sone()));
  check("Theta = x c dphi' - x", (grep { $th2->[$_] != $Theta->[$_] } 0..$N) == 0);
  check("Theta_0 = 0, Theta_1 = -x^3/c^2, Theta_2 = x^3 (1-x)^2/(2 c^4)",
        $Theta->[0] == 0 && $Theta->[1] == -$x**3/$c**2 && $Theta->[2] == $x**3*(1-$x)**2/(2*$c**4));
  # compare with the printed Theta(t), dphi(t) (polynomials of degree 7 in t): evaluate at 9 values of t
  for my $tv (map { R($_) } (1, 2, -1, "1/2", "-1/3", 3, "2/5", -2, "5/7")) {
    my ($pt, $pd) = (R(0), R(0)); for my $n (0..7) { $pt += $Theta->[$n]*$tv**$n; $pd += $dphi->[$n]*$tv**$n; }
    check("printed Theta(t) at t = $tv", meval($thetaL, $x, $tv) == $pt);
    check("printed dphi(t) at t = $tv", meval($dphiL, $x, $tv) == $pd);
  }
  my $uinv = sinv($u); my $w = szero(); $w->[$_] = $uinv->[$_+1] for 0..$N-1;     # 1/tau - 1/t (valid to t^(N-1))
  my $l = ssub(sscal($c, $w), sscal($x, $s));
  check("l_0 = c/2 - x/c", $l->[0] == $c/2 - $x/$c);
  # x^2 - c x L_2 = -x^3 tau s*^2 = Theta (L_2 = x s*), and (x - c L_1)(x - c L_2) = x^2 (multiplied by t; L_1 = c/t + l)
  my $th3 = ssub(sscal($x**2, sone()), sscal($c*$x*$x, $s));
  check("x^2 - c x L_2 = -x^3 tau s*^2", (grep { $th3->[$_] != $Theta->[$_] } 0..$N) == 0);
  my $tL1 = szero(); $tL1->[0] = $c; $tL1->[$_] = $l->[$_-1] for 1..$N;            # t L_1, valid to t^N
  my $prod = smul(ssub(sscal($x, [ map { $_ == 1 ? R(1) : R(0) } 0..$N ]), sscal($c, $tL1)), ssub(sscal($x, sone()), sscal($c*$x, $s)));
  check("(x - c L_1)(x - c L_2) = x^2 (times t), t^0..t^(N-1)", (grep { $prod->[$_] != ($_ == 1 ? $x**2 : 0) } 0..$N-1) == 0);
  # L_1 + L_2 = c/tau, L_1 L_2 = x/tau, times t: t L_1 + t L_2 = c/u, t L_1 L_2 = x/u
  my $sumL = sadd($tL1, [ map { $_ == 0 ? R(0) : $x*$s->[$_-1] } 0..$N ]);
  check("t (L_1 + L_2) = c/u", (grep { $sumL->[$_] != sscal($c, $uinv)->[$_] } 0..$N-1) == 0);
  check("t L_1 L_2 = x/u", (grep { smul($tL1, sscal($x,$s))->[$_] != sscal($x, $uinv)->[$_] } 0..$N-1) == 0);
  # alg and rho = tau/(t alg)
  my $ts = smul($tau, $s);
  my $p2 = sadd(sadd(sscal(-1, sinv(smul($s,$s))), sscal(-$x, sinv(smul(ssub(sone(),$s), ssub(sone(),$s))))),
                smul(smul($tau,$tau), sinv(smul(ssub(sone(),$ts), ssub(sone(),$ts)))));
  my $alg = smul(ssub(sone(), $ts), ssqrt1(sscal(-$x/$c**3, $p2)));
  my $rho = smul($u, sinv($alg));
  my @rhoD = map { meval($_, $x) } topsplit(firstlist($rhoL));
  check("rho = tau/(t alg) against the stored rho, t^0..t^5", (grep { $rho->[$_] != $rhoD[$_] } 0..$KD) == 0);
  my $rhoinv = sinv($rho);

  # ---------- dilaton data: b^(k)_j (k = 1..5), leaf coefficients a^(k)_j of g_k ----------
  my (@b, @a);
  $a[0] = [ R(0), $c ];                                         # g_0 = c lam_1
  for my $k (1..$KD) { my @v = map { meval($_, $x) } topsplit($dilL{$k}); $b[$k] = \@v;
    check("dilaton t^$k: no lam_0", $v[0] == 0);
    my $sum = R(0); $sum += $_ for @v; check("dilaton t^$k: Sum_j b_j = 0 (regular at z = 0)", $sum == 0);
    $a[$k] = [ map { -($_/$c)*$v[$_] } 0..$#v ];
    my $sa = R(0); $sa += $a[$k][$_]/$_ for 1..$#v; check("g_$k = O(psi^-2): Sum_j a_j/j = 0", $sa == 0);
    # top coefficient b^(k)_{k+1} = (-1)^(k+1) prod_{j=1}^k (k+1-jc)^2/(k!(k+1)! c^(2k-2))
    my $P = R(1); $P *= ($k+1-$_*$c)**2 for 1..$k;
    check("top dilaton coefficient t^$k", $v[$k+1] == R((-1)**($k+1))*$P/(fact($k)*fact($k+1)*$c**(2*$k-2)));
  }
  my $gval = sub { my ($k, $psi) = @_; my $r = R(0); for my $j (1..$#{$a[$k]}) { $r += $a[$k][$j]/($c - $j*$psi); } return $r; };
  # ===== residue-lemma control =====
  # sigma^(k) leaf coefficients by (6.4): comp0 lam_j: g_j^(k) = $a[k][j]; comp1 lam_j: (k+2-j) b_j^(k+1) + Sum_m l_m g_j^(k-m)  (b = q-leaf table, b^(0)_1 = -c^2)
  $b[0] = [ R(0), -$c**2 ];
  my %Y;   # $Y{n}[comp] = [const, lam_1, lam_2, ...] for Y = t sigma
  $Y{0} = [ [ R(0) ], [ $c ] ];
  for my $k (0..$KD-1) { my (@s0, @s1); $s0[0] = R(0); $s1[0] = R(0);
    for my $j (1..$k+2) { $s0[$j] = ($j <= $#{$a[$k]}) ? $a[$k][$j] : R(0);
      my $v = ($j <= $#{$b[$k+1]}) ? ($k+2-$j)*$b[$k+1][$j] : R(0);
      for my $m (0..$k) { $v += $l->[$m]*(($j <= $#{$a[$k-$m]}) ? $a[$k-$m][$j] : R(0)); } $s1[$j] = $v; }
    $Y{$k+1} = [ \@s0, \@s1 ]; }
  my $NY = $KD;   # Y^(0..NY) known
  my $pairc = sub { my ($A,$B) = @_; # coefficient matrix of <A(psi), B(psi')> in basis (1, lam_j) x (1, lam_l)
    my %M; my $add = sub { my ($u,$v,$f) = @_; for my $i (0..$#$u) { for my $ll (0..$#$v) { my $w = $f*($u->[$i]//R(0))*($v->[$ll]//R(0)); next if $w == 0; $M{"$i,$ll"} = ($M{"$i,$ll"}//R(0)) + $w; } } };
    $add->($A->[0],$B->[1],$x**2); $add->($A->[1],$B->[0],$x**2); $add->($A->[1],$B->[1],-$c*$x); return \%M; };
  my $kappa = -$x*$c**3;
  for my $n (0..$NY) { my %C;
    for my $n1 (0..$n) { my $M = $pairc->($Y{$n1}, $Y{$n-$n1}); $C{$_} = ($C{$_}//R(0)) - $M->{$_} for keys %$M; }
    my $Cg = sub { $C{"$_[0],$_[1]"} // R(0) }; my $maxi = 0; for (keys %C) { my ($i,$ll) = split /,/; $maxi = $i if $i > $maxi; $maxi = $ll if $ll > $maxi; }
    check("conifold residue: C^($n)_00 = -[t^$n] kappa", $Cg->(0,0) == ($n == 0 ? -$kappa : R(0)));
    for my $i (1..$maxi) { my $Sv = $c*$Cg->($i,0); $Sv += $i*$Cg->($i,$_)/($i+$_) for 1..$maxi; check("conifold residue: S_$i = 0 at t^$n (symplecticity)", $Sv == 0); }
    next if $n < 2; my $k = $n - 2; next if $k > 3;
    # prediction e^(k)_il = il/(c(i+l)) C^(k+2)_il versus the full loop matrices of the Mathematica tables (t^0..t^3, indices 0..k+1)
    my @rows = topsplit($loopFL{$k}); my $out = '';
    for my $i (0..$#rows) { my @col = topsplit($rows[$i]);
      for my $ll (0..$#col) { my $val = meval($col[$ll], $x); my $pred = ($i >= 1 && $ll >= 1) ? $i*$ll/($c*($i+$ll))*$Cg->($i,$ll) : R(0);
        check("conifold residue: e^($k)_($i,$ll) = il/(c(i+l)) C_il vs Mathematica loop table", $val == $pred); } }
    for my $i (1..$maxi) { for my $ll (1..$maxi) { next if ($i-1)+($ll-1) <= $k; check("conifold residue: support t^$k ($i,$ll)", $Cg->($i,$ll) == 0); } }
    if ($k <= 1) { my @e = map { my $i=$_; [ map { $i*$_/($c*($i+$_))*$Cg->($i,$_) } 1..$k+1 ] } 1..$k+1;
      print "  x = $x, [t^$k] rho^2 E11 by residue formula: ", join(' ; ', map { join(', ', @$_) } @e), "\n";
      if ($k == 0) { check("conifold residue: e^(0)_11 = -x^3/2", $e[0][0] == -$x**3/2); }
      if ($k == 1) { check("conifold residue: e^(1)_11 = -x^3(1-x)^2/(2c^2)", $e[0][0] == -$x**3*(1-$x)**2/(2*$c**2));
                     check("conifold residue: e^(1)_12 = 2x^3(1-x)^2/(3c^2)", $e[0][1] == 2*$x**3*(1-$x)**2/(3*$c**2));
                     check("conifold residue: e^(1)_22 = 0", $e[1][1] == 0); } } }
}
print "checks passed: $nok, failed: $nbad\n";
