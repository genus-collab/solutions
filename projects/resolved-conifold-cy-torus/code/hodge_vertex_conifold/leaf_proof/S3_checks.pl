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

my $dir = shift // '..';
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

  # ---------- closed form G(z; t) at rational z versus the dilaton tables (L1, z-form) ----------
  for my $zv (map { R($_) } ("3/11", "-2/19", 5, "-13/4")) {
    my $z = $zv;
    my $E = sexp(sscal(-1/$z, $dphi));
    my $F = szero(); my $tp = sone();
    for my $n (0..$N) { my $Fn = R(1); $Fn *= (1+$_*$z)**2 for 1..$n; my $den = fact($n)*$z**$n; $den *= ($c + $_*$z) for 2..$n+1;
      $F = sadd($F, sscal($Fn/$den, $tp)); $tp = smul($tp, $tau); }
    my $G = sscal($c/($c+$z), smul($u, smul($F, $E)));
    check("z G_0 = c - c^2/(c+z) at z = $z", $z*$G->[0] == $c - $c**2/($c+$z));
    for my $k (1..$KD) { my $pf = R(0); $pf += $b[$k][$_]/($c + $_*$z) for 1..$#{$b[$k]};
      check("[t^$k](z G) = Sum_j b_j/(c + j z) at z = $z (closed form vs dilaton table)", $z*$G->[$k] == $pf); }
    # G_1 in closed form: -(1-x)^2/(2 (c+z)(c+2z))
    check("G_1 = -(1-x)^2/(2(c+z)(c+2z)) at z = $z", $G->[1] == -(1-$x)**2/(2*($c+$z)*($c+2*$z)));
  }

  # ---------- leaf coefficients as functions of t (residues at z_j = -c/j) ----------
  # b_1(t) = -c^2 (tau/t) F(-c, tau) e^(dphi/c);  b_j(t) = j (c z_j/(c+z_j)) (tau/t) e^(j dphi/c) Sum_{n>=j-1} Res_{z_j}(F_n) tau^n
  for my $j (1..$KD+1) {
    my $zj = -$c/$j; my $Rs = szero(); my $tp = sone();
    for my $n (0..$N) { my $val = R(0);
      if ($j == 1 || $n >= $j-1) { my $num = R(1); $num *= (1 + $_*$zj)**2 for 1..$n; my $den = fact($n)*$zj**$n;
        for my $i (2..$n+1) { $den *= ($c + $i*$zj) unless $i == $j; } $val = ($j == 1) ? $num/$den : $num/($den*$j); }
      $Rs = sadd($Rs, sscal($val, $tp)); $tp = smul($tp, $tau); }
    my $pref = ($j == 1) ? -$c**2 : $j*$c*$zj/($c+$zj);
    my $bj = sscal($pref, smul($u, smul(sexp(sscal($j/$c, $dphi)), $Rs)));
    # hypergeometric form: b_j(t) = (-1)^j C_j t^(j-1) ((tau/t) e^(dphi/c))^j 2F1(jx/c, jx/c; j; tau), C_j = P_{j-1}/(j!(j-1)! c^(2j-4))
    my $P = R(1); $P *= ($j - $_*$c)**2 for 1..$j-1; my $Cj = $P/(fact($j)*fact($j-1)*$c**(2*$j-4));
    my $hyp = szero(); my $tq = sone(); my $coef = R(1);
    for my $m (0..$N) { $hyp = sadd($hyp, sscal($coef, $tq)); $tq = smul($tq, $tau); $coef *= ($j*$x/$c + $m)**2/(($j + $m)*($m + 1)); }
    my $om = smul($u, sexp(sscal(1/$c, $dphi))); my $omj = sone(); $omj = smul($omj, $om) for 1..$j;
    my $bh = sscal(R((-1)**$j)*$Cj, smul($omj, $hyp)); my $bh2 = szero(); $bh2->[$_ + $j - 1] = $bh->[$_] for 0..$N-$j+1;
    for my $k (0..$KD) { my $dat = ($k == 0) ? ($j == 1 ? -$c**2 : R(0)) : ($j <= $#{$b[$k]} ? $b[$k][$j] : R(0));
      check("leaf coefficient b_$j(t) by residues, t^$k", $bj->[$k] == $dat);
      check("leaf coefficient b_$j(t), hypergeometric form, t^$k", $bh2->[$k] == $dat); }
  }

  # ---------- loops: the integrated E11 flow ----------
  # e^(k)_{il} = (1/(k+2)) Sum_{k0>=1, k1+k2 = k+1-k0} Theta_k0 (k1+1)(k2+1) a^(k1)_i a^(k2)_l
  my @e;
  for my $k (0..$KD) { my %m;
    for my $k0 (1..$k+1) { for my $k1 (0..$k+1-$k0) { my $k2 = $k+1-$k0-$k1;
        for my $i (1..$k1+1) { for my $ll (1..$k2+1) {
            $m{"$i,$ll"} = ($m{"$i,$ll"} // R(0)) + $Theta->[$k0]*($k1+1)*($k2+1)*$a[$k1][$i]*$a[$k2][$ll]/($k+2); } } } }
    $e[$k] = \%m; }
  for my $k (0..$KD) {
    # support: nonzero only for (i-1)+(l-1) <= k
    check("loop t^$k: support (i-1)+(l-1) <= k", (grep { my ($i,$ll) = split /,/; ($i-1)+($ll-1) > $k && $e[$k]{$_} != 0 } keys %{$e[$k]}) == 0);
    check("loop t^$k: symmetric", (grep { my ($i,$ll) = split /,/; ($e[$k]{"$ll,$i"} // R(0)) != $e[$k]{$_} } keys %{$e[$k]}) == 0);
    # restricted data {i, l} -> e_il (i <= l)
    my %dat; for my $it (topsplit($loopRL{$k})) { $it =~ /^\{(\d+), (\d+)\} -> (.*)$/ or die "bad item $it\n"; $dat{"$1,$2"} = meval($3, $x); }
    for my $key (keys %{$e[$k]}) { my ($i,$ll) = split /,/, $key; next if $i > $ll;
      check("loop t^$k e_$i$ll (formula vs restricted table)", $e[$k]{$key} == ($dat{$key} // R(0))); }
    for my $key (keys %dat) { check("loop t^$k: table entry $key covered", exists $e[$k]{$key}); }
    # corners (i-1)+(l-1) = k: (-1)^(k+1) x^3 (k1+1)(k2+1) P_k1 P_k2/((k+2) k1!^2 k2!^2 c^(2k))
    for my $k1 (0..$k) { my $k2 = $k - $k1; my ($P1, $P2) = (R(1), R(1));
      $P1 *= ($k1+1-$_*$c)**2 for 1..$k1; $P2 *= ($k2+1-$_*$c)**2 for 1..$k2;
      my $pred = R((-1)**($k+1))*$x**3*($k1+1)*($k2+1)*$P1*$P2/(($k+2)*fact($k1)**2*fact($k2)**2*$c**(2*$k));
      check("loop t^$k corner (".($k1+1).",".($k2+1).")", ($e[$k]{($k1+1).",".($k2+1)} // R(0)) == $pred); }
    if ($k <= 3) {   # full matrices (lam_0..lam_{k+1}) including the zeros
      my @rows = topsplit($loopFL{$k});
      for my $i (0..$#rows) { my @col = topsplit($rows[$i]);
        for my $ll (0..$#col) { my $val = meval($col[$ll], $x); my $pred = ($i >= 1 && $ll >= 1) ? ($e[$k]{"$i,$ll"} // R(0)) : R(0);
          check("full loop matrix t^$k ($i,$ll)", $val == $pred); } } }
  }
  check("rank one at t^0: e^(0) = -(x^3/2) lam_1 lam_1'", scalar(keys %{$e[0]}) == 1 && $e[0]{"1,1"} == -$x**3/2);

  # ---------- legs j = 1: [t^k] rho f_{1,b} = Sum_m rhoinv_m Sum_l e^(k-m)_{il} l^b/c^(b+1) ----------
  for my $bb (0..7) { for my $k (0..$KD) { my $key = "1,$bb,$k"; next unless exists $legL{$key};
      my @v = map { meval($_, $x) } topsplit($legL{$key});
      my @pred = map { R(0) } 0..$k+1;
      for my $m (0..$k) { for my $kk (keys %{$e[$k-$m]}) { my ($i,$ll) = split /,/, $kk;
          $pred[$i] += $rhoinv->[$m]*$e[$k-$m]{$kk}*R($ll)**$bb/$c**($bb+1); } }
      check("leg j = 1, b = $bb, t^$k", (grep { $v[$_] != $pred[$_] } 0..$k+1) == 0);
      if ($k == 0) { check("leg t^0 closed form, b = $bb", $v[1] == -$x**3/(2*$c**($bb+1)) && $v[0] == 0); }
      if ($k == 1) { check("leg t^1 closed form, b = $bb", $v[1] == (1-$x)**2*$x**3*(R(2)**($bb+3) - 9)/(12*$c**($bb+3))
                           && $v[2] == (1-$x)**2*$x**3*8/(12*$c**($bb+3)) && $v[0] == 0); } } }

  # ---------- sigma: symplectic check and loops directly from sigma ----------
  for my $pp (["1/13", "1/23"], ["-2/17", "-5/29"], ["3/19", "4/31"]) {
    my ($psi, $psp) = (R($pp->[0]), R($pp->[1]));
    my $tsig = sub { my $p = shift; my ($s0, $s1) = (szero(), szero());
      for my $m (0..$KD) { $s0->[$m+1] = $gval->($m, $p) if $m+1 <= $N;
        my $v = ($c - ($m+1)*$p)*$gval->($m, $p); $v += $l->[$_]*$gval->($m-1-$_, $p) for 0..$m-1; $s1->[$m] = $v; }
      return ($s0, $s1); };
    my $pair = sub { my ($a0,$a1,$b0,$b1) = @_; return ssub(sscal($x**2, sadd(smul($a0,$b1), smul($a1,$b0))), sscal($c*$x, smul($a1,$b1))); };
    my ($A0, $A1) = $tsig->($psi); my ($B0, $B1) = $tsig->(-$psi); my ($C0, $C1) = $tsig->($psp);
    check("[t^0] t sigma = (0, c) at psi = $psi", $A0->[0] == 0 && $A1->[0] == $c);
    check("[t^0] sigma_1 = x lam_1 at psi = $psi", $A1->[1] == $x/($c - $psi));
    check("[t^1] sigma_1 = (x/(2c^2))[(1-x+x^2) lam_1 - 2(1-x)^2 lam_2] at psi = $psi",
          $A1->[2] == ($x/(2*$c**2))*((1-$x+$x**2)/($c-$psi) - 2*(1-$x)**2/($c-2*$psi)));
    check("[t^1] g = ((1-x)^2/(2c))(lam_1 - 2 lam_2) at psi = $psi", $A0->[2] == ((1-$x)**2/(2*$c))*(1/($c-$psi) - 2/($c-2*$psi)));
    my $S = $pair->($A0, $A1, $B0, $B1);
    check("symplecticity <t sigma(psi), t sigma(-psi)> = -x c^3 + O(t^6) at psi = $psi",
          $S->[0] == -$x*$c**3 && (grep { $S->[$_] != 0 } 1..$KD) == 0);
    # pi(sigma) = <sigma, (1, L_2)> = x^3 tau s*^2 (psi/t) d/dt(t g):  check t*pi(t sigma)/t ... at orders t^0..t^4
    my $L2 = sscal($x, $s);
    my $pis = sadd(smul(sscal($x**2, $L2), $A0), smul(ssub(sscal($x**2, sone()), sscal($c*$x, $L2)), $A1));   # <t sigma, (1, L_2)>
    my $dtg = szero(); $dtg->[$_] = ($_+1)*$gval->($_, $psi) for 0..$KD;                                     # d/dt (t g)
    my $pred = sscal($psi, smul(sscal(-1, $Theta), $dtg));                                                    # t pi(sigma) = -Theta psi d/dt(t g)
    check("pi(sigma) = -(x^2 - c x L_2)(psi/t) d/dt(t g), t^0..t^5 at psi = $psi", (grep { $pis->[$_] != $pred->[$_] } 0..$KD) == 0);
    # loops: t^2 rho^2 E11 = [-x c^3 - <t sigma(psi), t sigma(psi')>]/(psi + psi')
    my $Np = $pair->($A0, $A1, $C0, $C1);
    my $Phi = sscal(-1/($psi+$psp), sadd(sscal($x*$c**3, sone()), $Np));
    check("t^2 rho^2 E11 = O(t^2) at ($psi, $psp)", $Phi->[0] == 0 && $Phi->[1] == 0);
    for my $k (0..$KD-2) { my $v = R(0);
      for my $kk (keys %{$e[$k]}) { my ($i,$ll) = split /,/, $kk; $v += $e[$k]{$kk}/(($c - $i*$psi)*($c - $ll*$psp)); }
      check("[t^$k] rho^2 E11 from sigma = integrated flow at ($psi, $psp)", $Phi->[$k+2] == $v); }
  }
}
print "checks passed: $nok, failed: $nbad\n";
