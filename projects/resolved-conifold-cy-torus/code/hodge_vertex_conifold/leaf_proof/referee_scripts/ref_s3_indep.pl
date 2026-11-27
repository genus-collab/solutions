#!/usr/bin/perl
# ref_s3_indep.pl -- referee's independent exact checks for leaf_proof/S3_leaves_loops.md (Perl + Math::BigRat).
# Independent of S3_checks.pl: the frame-B singular rows are built from the closed form (S1: G = rho r_{0,1}) as exact
# double series in (t, z) INCLUDING the negative z-powers (so regularity at z = 0 is tested, not assumed); rho from
# Delta_1 = x^2 W2(U1) (V1); the H-row from the code definition (rmatrixCore); E_11 from the code's edge formula (evalQ:
# EdN = Sum_k (-1)^k pre(m+k+1, n-k)) -- NOT from (11)/(12) -- and then compared with (i) the integrated flow (16),
# (ii) the stored Mathematica loop / leg / dilaton tables, at x values not used by the author.
# Usage: perl ref_s3_indep.pl <project dir>
use strict; use warnings; use Math::BigRat;

my $dir = shift // '.';
my $NT = 8; my $MZ = 9; my $ZLO = -$NT; my $ZHI = $MZ + $NT; my $NZ = $ZHI - $ZLO + 1;
sub Bq { Math::BigRat->new("$_[0]") }
my ($nok, $nbad) = (0, 0); my %cat;
sub check { my ($cat, $ok, $msg) = @_; $cat{$cat}++; if ($ok) { $nok++ } else { $nbad++; print "FAIL [$cat] $msg\n"; } }

# ---------------- t-series 0..NT
sub tz { [ map { Bq(0) } 0..$NT ] }
sub tc { my $s = tz(); $s->[0] = Bq($_[0]); $s }
sub tadd { my ($a,$b)=@_; [ map { $a->[$_] + $b->[$_] } 0..$NT ] }
sub tsub { my ($a,$b)=@_; [ map { $a->[$_] - $b->[$_] } 0..$NT ] }
sub tsc { my ($k,$a)=@_; [ map { $a->[$_] * $k } 0..$NT ] }
sub tmul { my ($a,$b)=@_; my $r = tz(); for my $i (0..$NT) { next if $a->[$i]->is_zero; for my $j (0..$NT-$i) { $r->[$i+$j] += $a->[$i]*$b->[$j]; } } $r }
sub tinv { my $a=shift; die "tinv\n" if $a->[0]->is_zero; my $r=tz(); $r->[0] = 1/$a->[0];
  for my $n (1..$NT) { my $s=Bq(0); $s += $a->[$_]*$r->[$n-$_] for 1..$n; $r->[$n] = -$s/$a->[0]; } $r }
sub tder { my $a=shift; my $r=tz(); $r->[$_-1] = $a->[$_]*$_ for 1..$NT; $r }          # valid to t^(NT-1)
sub ttd { my $a=shift; [ map { $a->[$_]*$_ } 0..$NT ] }                                   # t d/dt
sub tint { my $a=shift; my $r=tz(); $r->[$_] = $a->[$_-1]/$_ for 1..$NT; $r }
sub tshift { my ($a,$s)=@_; my $r=tz(); for my $i (0..$NT) { my $j=$i+$s; $r->[$j] = $a->[$i] if $j>=0 && $j<=$NT; } $r }  # times t^s (s may be <0: drops)
sub tlog1 { my $a=shift; die "tlog1\n" unless $a->[0] == 1; tint(tmul(tder($a), tinv($a))) }
sub texp0 { my $a=shift; die "texp0\n" unless $a->[0]->is_zero; my $r=tc(1);
  for my $n (1..$NT) { my $s=Bq(0); $s += $a->[$_]*$r->[$n-$_]*$_ for 1..$n; $r->[$n] = $s/$n; } $r }
sub tpow1 { my ($a,$p)=@_; texp0(tsc($p, tlog1($a))) }                                      # a^p, a(0) = 1
sub teq { my ($a,$b,$upto)=@_; for my $i (0..$upto) { return 0 unless $a->[$i] == $b->[$i]; } 1 }
sub fact { my $n=shift; my $r=Bq(1); $r *= $_ for 1..$n; $r }

# ---------------- z power series 0..LZ-1 (for the Taylor parts of F_n)
my $LZ = $ZHI + $NT + 2;
sub zser_lin_inv { my ($a,$b)=@_; [ map { (Bq(-$b)/$a)**$_ / $a } 0..$LZ-1 ] }               # 1/(a + b z)
sub zmul { my ($p,$q)=@_; my @r = map { Bq(0) } 0..$LZ-1; for my $i (0..$LZ-1) { next if $p->[$i]->is_zero; for my $j (0..$LZ-1-$i) { $r[$i+$j] += $p->[$i]*$q->[$j]; } } \@r }

# ---------------- double series: d->[k][j - ZLO], k = 0..NT, j = ZLO..ZHI
sub dz { [ map { [ map { Bq(0) } 1..$NZ ] } 0..$NT ] }
sub dmul { my ($a,$b,$zmax)=@_; my $r = dz();
  for my $k1 (0..$NT) { for my $i1 (0..$NZ-1) { my $v1 = $a->[$k1][$i1]; next if $v1->is_zero;
    for my $k2 (0..$NT-$k1) { for my $i2 (0..$NZ-1) { my $v2 = $b->[$k2][$i2]; next if $v2->is_zero;
      my $j = $i1 + $i2 + 2*$ZLO; next if $j > $zmax; die "below ZLO\n" if $j < $ZLO; $r->[$k1+$k2][$j-$ZLO] += $v1*$v2; } } } }
  $r }

# ---------------- parsing of Mathematica InputForm
sub mconv { my $s = shift; $s =~ s/\^(\d+)/**$1/g; $s =~ s/\bx\b/\$x/g; $s =~ s/(?<!\*\*)(?<![\d.])(\d+)/Bq($1)/g; $s }
sub meval { my ($s,$xv) = @_; my $x = $xv; my $v = eval mconv($s); die "eval failed ($@): $s\n" if $@; Bq($v) }
sub topsplit { my $s = shift; $s =~ s/^\s*\{//; $s =~ s/\}\s*$//; my @out; my $d = 0; my $cur = '';
  for my $ch (split //, $s) { $d++ if $ch eq '(' || $ch eq '{'; $d-- if $ch eq ')' || $ch eq '}';
    if ($ch eq ',' && $d == 0) { push @out, $cur; $cur = ''; } else { $cur .= $ch; } } push @out, $cur if $cur =~ /\S/;
  s/^\s+|\s+$//g for @out; @out }
sub firstlist { my $s = shift; my $i = index($s, '{'); my $d = 0; for my $k ($i..length($s)-1) { my $ch = substr($s,$k,1);
    $d++ if $ch eq '{'; $d-- if $ch eq '}'; return substr($s, $i, $k-$i+1) if $d == 0; } die "unbalanced\n"; }
sub slurp { my $f = shift; open(my $fh, '<', $f) or die "cannot open $f\n"; my @l = <$fh>; close $fh; chomp @l; @l }
my @kt5 = slurp("$dir/hvc_leaf_structure_out_kt5.txt");
my ($rhoL) = grep { /^rho = / } @kt5;
my (%dilL, %legL, %loopL);
for (@kt5) {
  $dilL{$1} = firstlist($2) if /^\(i\) rho q_t - q_0, t\^(\d+): coefficients of lam_0..lam_\d+ = (.*)$/;
  $legL{"$1,$2,$3"} = firstlist($4) if /^\(ii\) leg j = (\d), b = (\d+), t\^(\d+): lam_0..lam_\d+ = (\{.*)$/;
  $loopL{$1} = firstlist($2) if /^\(iii'\) rho\^2 E11, t\^(\d+), restricted ansatz \(no lam_0, D <= k\): (.*)$/;
}

for my $xs ("3/7", "-2/5", "5/2", "4/3") {
  my $x = Bq($xs); my $c = 1 + $x;
  print "---- x = $x\n";
  # tau, u = tau/t, t/tau, Q = e^t
  my $tau = tz(); $tau->[$_] = Bq((-1)**($_+1))/fact($_) for 1..$NT;
  my $u = tz(); $u->[$_] = Bq((-1)**$_)/fact($_+1) for 0..$NT;
  my $ui = tinv($u);
  my $Q = tz(); $Q->[$_] = 1/fact($_) for 0..$NT;
  # s*: root of x tau s^2 - c s + 1 with s(0) = 1/c
  my $s = tc(1/$c); for (0..$NT+2) { $s = tsc(1/$c, tadd(tc(1), tsc($x, tmul($tau, tmul($s,$s))))); }
  my $ts = tmul($tau, $s);
  # dphi from its DEFINITION phi_tau(s*) - phi_0(1/c) (logs of series with constant term 1)
  my $dphi = tadd(tadd(tlog1(tsc($c, $s)), tsc($x, tlog1(tsc($c/$x, tsub(tc(1), $s))))), tsc(-1, tlog1(tsub(tc(1), $ts))));
  check("dphi' = 1 - x s* (S1 1.2(c), S3 (14))", teq(tder($dphi), tsub(tc(1), tsc($x,$s)), $NT-1), "");
  # U1, U2, W2, Delta_1, Delta_2 from the code's superpotential (V1)
  my $U1 = tsub(tc(1), $ts); my $U2 = tsub(tc(($x-1)/$x), $U1);
  # t^2 Delta_1 = -x^2 (t/tau)^2 (1 - tau s*) [x (1 - tau)/(1 - s*)^2 + 1/s*^2]   (W2 = (U d/dU)^2 W at U1)
  my $om = tsub(tc(1), $s);
  my $br = tadd(tsc($x, tmul(tsub(tc(1),$tau), tinv(tmul($om,$om)))), tinv(tmul($s,$s)));
  my $t2D1 = tsc(-$x*$x, tmul(tmul($ui,$ui), tmul($U1, $br)));
  check("t^2 Delta_1 (0) = -x c^3", $t2D1->[0] == -$x*$c**3, "");
  # cross-check the W2 formula against a direct evaluation of -xQU/(1-QU)^2 - U/(1-U)^2 times tau^2 (both regular after x tau^2)
  my $QU = tmul($Q, $U1); my $omQU = tsub(tc(1), $QU);                                   # 1 - Q U1 = O(t): divide by t
  my $omQUt = tshift($omQU, -1); my $omUt = tshift(tsub(tc(1), $U1), -1);                  # valid to t^(NT-1)
  my $W2t2 = tadd(tsc(-$x, tmul($QU, tinv(tmul($omQUt,$omQUt)))), tsc(-1, tmul($U1, tinv(tmul($omUt,$omUt)))));   # t^2 W2
  check("t^2 W2(U1) direct = closed formula", teq(tsc($x*$x, $W2t2), $t2D1, $NT-2), "");
  my $rho = tpow1(tsc(1/(-$x*$c**3), $t2D1), Bq("-1/2"));
  my $rhoi = tinv($rho);
  my @rhoD = map { meval($_, $x) } topsplit(firstlist($rhoL));
  check("rho vs stored table t^0..t^5", (grep { $rho->[$_] != $rhoD[$_] } 0..5) == 0, "rho");
  my $QU2 = tmul($Q, $U2); my $D2 = tsc($x*$x, tadd(tsc(-$x, tmul($QU2, tinv(tmul(tsub(tc(1),$QU2), tsub(tc(1),$QU2))))), tsc(-1, tmul($U2, tinv(tmul(tsub(tc(1),$U2), tsub(tc(1),$U2)))))));
  check("Delta_2(0) = x^3/c", $D2->[0] == $x**3/$c, "");
  # eigenvalues from LofU: L = (x U + 1)/tau;  t L_1 = (t/tau)(x U1 + 1),  L_2 = (x U2 + 1)/tau with x U2 + 1 = O(t)
  my $tL1 = tmul($ui, tadd(tsc($x,$U1), tc(1)));
  my $L2 = tmul($ui, tshift(tadd(tsc($x,$U2), tc(1)), -1));                               # valid to t^(NT-1)
  check("L_2 = x s* (from LofU)", teq($L2, tsc($x,$s), $NT-1), "");
  # (9): Delta_2 = (L_2 - L_1)(x^2 - c x L_2),  t^2 Delta_1 = t(L_1 - L_2) t(x^2 - c x L_1)
  my $tdL = tsub(tshift($L2,1), $tL1);                                                    # t (L_2 - L_1), valid to t^NT? L2 to NT-1 -> shifted to NT
  my $B = tsub(tc($x*$x), tsc($c*$x, $L2));                                                # x^2 - c x L_2 = O(t)
  my $Bt = tshift($B, -1);                                                                 # (x^2 - c x L_2)/t
  check("(9) Delta_2 = (L_2 - L_1)(x^2 - c x L_2)", teq(tmul($tdL, $Bt), $D2, $NT-3), "");
  my $tA1 = tsub(tc(0), $tdL); my $C1 = tsub(tshift(tc($x*$x),1), tsc($c*$x, $tL1));      # t (x^2 - c x L_1)
  check("(9) t^2 Delta_1 = t(L_1 - L_2) t(x^2 - c x L_1)", teq(tmul($tA1, $C1), $t2D1, $NT-3), "");
  # Theta: V8 definition (L_2 - L_1)(c x L_2 - x^2)^2/Delta_2 versus -x^3 tau s*^2 and x c dphi' - x
  my $Theta = tmul(tmul($tdL, tmul($B, $Bt)), tinv($D2));                                   # (t(L2-L1)) (B) (B/t) / Delta_2
  my $Th2 = tsc(-$x**3, tmul($tau, tmul($s,$s)));
  my $Th3 = tsub(tsc($x*$c, tder($dphi)), tc($x));
  check("Theta (V8 def) = -x^3 tau s*^2", teq($Theta, $Th2, $NT-3), "");
  check("Theta = x c dphi' - x", teq($Th2, $Th3, $NT-1), "");
  check("Theta = x^2 - c x L_2", teq($Th2, $B, $NT-1), "");

  # ---------------- closed form G(z; t) as a double series, with the negative z-powers kept
  # F = Sum_n tau^n z^-n T_n(z),  T_n = prod_{i=1}^n (1 + i z)^2 / (n! prod_{i=2}^{n+1} (c + i z))
  my $F = dz(); my $taun = tc(1);
  for my $n (0..$NT) {
    my $T = [ map { Bq(0) } 0..$LZ-1 ]; $T->[0] = 1/fact($n);
    for my $i (1..$n) { my $l = [ map { Bq(0) } 0..$LZ-1 ]; $l->[0] = Bq(1); $l->[1] = Bq($i); $T = zmul($T, zmul($l, $l)); }
    for my $i (2..$n+1) { $T = zmul($T, zser_lin_inv($c, Bq($i))); }
    for my $k (0..$NT) { next if $taun->[$k]->is_zero; for my $m (0..$LZ-1) { my $j = $m - $n; next if $j > $ZHI; $F->[$k][$j-$ZLO] += $taun->[$k]*$T->[$m]; } }
    $taun = tmul($taun, $tau);
  }
  # Exp[-dphi/z] = Sum_l (-1)^l dphi^l z^-l / l!
  my $E = dz(); my $dpl = tc(1);
  for my $l (0..$NT) { for my $k (0..$NT) { $E->[$k][-$l-$ZLO] += Bq((-1)**$l)*$dpl->[$k]/fact($l); } $dpl = tmul($dpl, $dphi); }
  my $H = dmul($F, $E, $MZ);
  for my $k (0..$NT) { check("regularity at z = 0 of [t^k](F Exp[-dphi/z]) (S1 Cor 7.1), independent", (grep { !$H->[$k][$_-$ZLO]->is_zero } $ZLO..-1) == 0, "k=$k"); }
  # times c/(c+z) and u = tau/t
  my $cz = dz(); $cz->[0][$_-$ZLO] = Bq((-1)**$_)/$c**$_ for 0..$ZHI;
  my $ud = dz(); $ud->[$_][-$ZLO] = $u->[$_] for 0..$NT;
  my $G = dmul(dmul($H, $cz, $MZ), $ud, $MZ);
  for my $k (0..$NT) { check("regularity at z = 0 of [t^k] G (S1 Cor 7.1), independent", (grep { !$G->[$k][$_-$ZLO]->is_zero } $ZLO..-1) == 0, "k=$k"); }
  my @Gk; for my $k (0..$NT) { $Gk[$k] = [ map { $G->[$k][$_-$ZLO] } 0..$MZ ]; }
  # G_0 = c/(c+z)
  check("F2: G_0 = c/(c+z)", (grep { $Gk[0][$_] != Bq((-1)**$_)/$c**$_ } 0..$MZ) == 0, "");

  # ---------------- dilaton tables: [psi^(m+1)] rho q_t = -(-1)^m G_k[m] versus Sum_j b_j j^(m+1)/c^(m+2)
  for my $k (1..5) { my @b = map { meval($_, $x) } topsplit($dilL{$k});
    my $ok = 1; for my $m (0..$MZ) { my $pred = Bq(0); $pred += $b[$_]*Bq($_)**($m+1)/$c**($m+2) for 1..$#b; $pred += $b[0]/$c if 0;    # b[0] = 0 (lam_0)
      $ok = 0 unless $pred == -Bq((-1)**$m)*$Gk[$k][$m]; }
    check("dilaton table t^k vs closed form, psi^1..psi^10", $ok && $b[0] == 0, "k=$k");
    my $sb = Bq(0); $sb += $_ for @b; check("Sum_j b_j = 0", $sb == 0, "k=$k");
    my $P = Bq(1); $P *= ($k+1-$_*$c)**2 for 1..$k;
    check("top coefficient (5)", $b[$k+1] == Bq((-1)**($k+1))*$P/(fact($k)*fact($k+1)*$c**(2*$k-2)), "k=$k"); }
  # hypergeometric form of the leaf coefficients (1.6): b_j(t) = (-1)^j C_j t^(j-1) ((tau/t) e^(dphi/c))^j 2F1(jx/c, jx/c; j; tau)
  { my @btab; $btab[0] = [ Bq(0), -$c**2 ]; for my $k (1..5) { $btab[$k] = [ map { meval($_, $x) } topsplit($dilL{$k}) ]; }
    my $om = tmul($u, texp0(tsc(1/$c, $dphi)));
    for my $j (1..6) { my $P = Bq(1); $P *= ($j - $_*$c)**2 for 1..$j-1; my $Cj = $P/(fact($j)*fact($j-1)*$c**(2*$j-4));
      my $hyp = tz(); my $tq = tc(1); my $co = Bq(1);
      for my $m (0..$NT) { $hyp = tadd($hyp, tsc($co, $tq)); $tq = tmul($tq, $tau); $co *= ($j*$x/$c + $m)**2/(($j + $m)*($m + 1)); }
      my $omj = tc(1); $omj = tmul($omj, $om) for 1..$j;
      my $bj = tshift(tsc(Bq((-1)**$j)*$Cj, tmul($omj, $hyp)), $j-1);
      for my $k (0..5) { my $dat = ($j <= $#{$btab[$k]}) ? $btab[$k][$j] : Bq(0);
        check("hypergeometric b_j(t) (1.6) vs dilaton table", $bj->[$k] == $dat, "j=$j k=$k"); } } }

  # ---------------- Ytilde = t rho v_1(z): (t G, (t L_1 + z) G + t z dG/dt), z-coefficients as t-series
  my (@Y0, @Y1);
  for my $m (0..$MZ) { my $g = [ map { $Gk[$_][$m] } 0..$NT ];
    $Y0[$m] = tshift($g, 1);
    my $v = tmul($tL1, $g);
    if ($m >= 1) { my $gp = [ map { $Gk[$_][$m-1] } 0..$NT ]; $v = tadd($v, tadd($gp, ttd($gp))); }
    $Y1[$m] = $v; }
  my $pair = sub { my ($a0,$a1,$b0,$b1) = @_; tsub(tsc($x*$x, tadd(tmul($a0,$b1), tmul($a1,$b0))), tsc($c*$x, tmul($a1,$b1))) };
  # Prop 4(a): [t^0] Ytilde = (0, c) (all z-orders)
  check("Prop 4(a): [t^0] t sigma = (0, c)", (grep { !$Y0[$_]->[0]->is_zero || $Y1[$_]->[0] != ($_ == 0 ? $c : 0) } 0..$MZ) == 0, "");
  # sigma^(0)_1 = x lam_1, sigma^(1) explicit (psi-form: [psi^m] = (-1)^m [z^m]; lam_j(psi) -> j^m/c^(m+1))
  my $lam = sub { my ($j,$m) = @_; Bq($j)**$m/$c**($m+1) };
  check("sigma^(0)_1 = x lam_1", (grep { Bq((-1)**$_)*$Y1[$_]->[1] != $x*$lam->(1,$_) } 0..$MZ) == 0, "");
  check("sigma^(1)_1 explicit", (grep { Bq((-1)**$_)*$Y1[$_]->[2] != ($x/(2*$c**2))*((1-$x+$x*$x)*$lam->(1,$_) - 2*(1-$x)**2*$lam->(2,$_)) } 0..$MZ) == 0, "");
  check("sigma^(1)_0 = g_1 explicit", (grep { Bq((-1)**$_)*$Y0[$_]->[2] != ((1-$x)**2/(2*$c))*($lam->(1,$_) - 2*$lam->(2,$_)) } 0..$MZ) == 0, "");
  check("l_0 = c/2 - x/c, l_1 = c/12 - x^2/c^3", $tL1->[1] == $c/2 - $x/$c && $tL1->[2] == $c/12 - $x*$x/$c**3, "");
  # symplecticity (I6): Sum_q (-1)^q <Yt[n-q], Yt[q]> = delta_n0 (-x c^3)
  for my $n (0..$MZ) { my $S = tz(); for my $q (0..$n) { $S = tadd($S, tsc(Bq((-1)**$q), $pair->($Y0[$n-$q], $Y1[$n-$q], $Y0[$q], $Y1[$q]))); }
    check("K_1 = 1 on the closed form, w^n", teq($S, ($n == 0 ? tc(-$x*$c**3) : tz()), $NT), "n=$n"); }
  # QDE row 2 (I4): d/dt Yt_1[m-1] = -x (t/tau) G[m] + L_2 Yt_1[m]   (m >= 1)
  for my $m (1..$MZ) { my $g = [ map { $Gk[$_][$m] } 0..$NT ];
    my $lhs = tder($Y1[$m-1]); my $rhs = tadd(tsc(-$x, tmul($ui, $g)), tmul($L2, $Y1[$m]));
    check("QDE row 2 on the closed form (I4)", teq($lhs, $rhs, $NT-2), "m=$m"); }
  # QDE row 1 sanity (definition): d/dt Yt_0[m-1] = Yt_1[m] - L_1 Yt_0[m]  (t L_1 Yt_0/t)
  for my $m (1..$MZ) { my $lhs = tder($Y0[$m-1]); my $rhs = tsub($Y1[$m], tmul($tL1, tshift($Y0[$m], -1)));
    check("QDE row 1 (definition)", teq($lhs, $rhs, $NT-2), "m=$m"); }

  # ---------------- E_11 from the CODE's edge formula: Phi(m,n) = t^2 rho^2 EdN(m,n) = Sum_k (-1)^k pre(m+k+1, n-k)
  my %pre; for my $p (0..$MZ) { for my $q (0..$MZ) { $pre{"$p,$q"} = tsc(-Bq((-1)**($p+$q)), $pair->($Y0[$p], $Y1[$p], $Y0[$q], $Y1[$q])); } }
  my %Phi;
  for my $m (0..$MZ-1) { for my $n (0..$MZ-1-$m) { my $v = tz(); for my $k (0..$n) { $v = tadd($v, tsc(Bq((-1)**$k), $pre{($m+$k+1).",".($n-$k)})); } $Phi{"$m,$n"} = $v; } }
  check("F3: Phi = t^2 rho^2 E_11 = O(t^2) (code formula)", (grep { !$Phi{$_}->[0]->is_zero || !$Phi{$_}->[1]->is_zero } keys %Phi) == 0, "");
  check("E_11(a,b) = E_11(b,a) (code formula)", (grep { my ($m,$n) = split /,/; !teq($Phi{$_}, $Phi{"$n,$m"}, $NT) } keys %Phi) == 0, "");
  # (16): [t^(k+2)] Phi = (1/(k+2)) Sum_{k0>=1, k0+k1+k2=k+1} Theta_k0 (k1+1)(k2+1) g_k1[m] g_k2[n],  g_k[m] = (-1)^m G_k[m]
  for my $key (keys %Phi) { my ($m,$n) = split /,/, $key; my $ok = 1;
    for my $k (0..$NT-3) { my $pred = Bq(0);
      for my $k0 (1..$k+1) { for my $k1 (0..$k+1-$k0) { my $k2 = $k+1-$k0-$k1;
          $pred += $Theta->[$k0]*($k1+1)*($k2+1)*Bq((-1)**($m+$n))*$Gk[$k1][$m]*$Gk[$k2][$n]; } }
      $pred /= ($k+2); $ok = 0 unless $Phi{$key}->[$k+2] == $pred; }
    check("Theorem 7 / (16) vs code edge formula, t^0..t^(NT-2)", $ok, $key); }
  # (16) uses Theta_k0 up to k0 = NT-1: Theta valid to t^(NT-3) only -> recheck with Th2 (valid to NT-1)
  for my $key (keys %Phi) { my ($m,$n) = split /,/, $key; my $ok = 1;
    for my $k (0..$NT-2) { my $pred = Bq(0);
      for my $k0 (1..$k+1) { for my $k1 (0..$k+1-$k0) { my $k2 = $k+1-$k0-$k1;
          $pred += $Th2->[$k0]*($k1+1)*($k2+1)*Bq((-1)**($m+$n))*$Gk[$k1][$m]*$Gk[$k2][$n]; } }
      $pred /= ($k+2); $ok = 0 unless $Phi{$key}->[$k+2] == $pred; }
    check("(16) with Theta = -x^3 tau s*^2", $ok, $key); }
  # loop tables: [t^k] rho^2 E_11 (m,n) = Sum e_il i^m l^n / c^(m+n+2)   (restricted symmetric ansatz, i <= l stored)
  for my $k (0..5) { my %e; for my $it (topsplit($loopL{$k})) { $it =~ /^\{(\d+), (\d+)\} -> (.*)$/ or die "bad $it\n"; $e{"$1,$2"} = meval($3, $x); }
    for my $key (keys %Phi) { my ($m,$n) = split /,/, $key; my $pred = Bq(0);
      for my $il (keys %e) { my ($i,$l) = split /,/, $il; my $v = $e{$il}*Bq($i)**$m*Bq($l)**$n/$c**($m+$n+2);
        $v += $e{$il}*Bq($l)**$m*Bq($i)**$n/$c**($m+$n+2) if $i != $l; $pred += $v; }
      check("loop table t^k vs code edge formula", $Phi{$key}->[$k+2] == $pred, "k=$k ($key)"); }
    # corner (17) and support of the table
    for my $k1 (0..$k) { my $k2 = $k - $k1; my ($i, $l) = ($k1+1, $k2+1); my $key = $i <= $l ? "$i,$l" : "$l,$i";
      my ($P1, $P2) = (Bq(1), Bq(1)); $P1 *= ($k1+1-$_*$c)**2 for 1..$k1; $P2 *= ($k2+1-$_*$c)**2 for 1..$k2;
      my $pred = Bq((-1)**($k+1))*$x**3*($k1+1)*($k2+1)*$P1*$P2/(($k+2)*fact($k1)**2*fact($k2)**2*$c**(2*$k));
      check("corner (17) vs loop table", ($e{$key} // Bq(0)) == $pred, "k=$k ($i,$l)"); }
    check("support (i-1)+(l-1) <= k of loop table", (grep { my ($i,$l) = split /,/; $i + $l - 2 > $k } keys %e) == 0, "k=$k"); }
  # leg tables j = 1: [t^k] rho E_11(a, b) = [t^(k+2)] (rho^-1 Phi(a,b)) vs Sum_i coef_i i^a / c^(a+1)
  for my $kb (keys %legL) { my ($j, $bb, $k) = split /,/, $kb; next unless $j == 1; next if $legL{$kb} =~ /Missing/;
    my @co = map { meval($_, $x) } topsplit($legL{$kb});
    for my $a (0..$MZ-1-$bb) { my $val = tmul($rhoi, $Phi{"$a,$bb"})->[$k+2];
      my $pred = Bq(0); $pred += $co[$_]*Bq($_)**$a/$c**($a+1) for 1..$#co; $pred += $co[0]/$c if $a == 0;
      check("leg table j=1 vs code edge formula", $val == $pred, "b=$bb k=$k a=$a"); } }
}
print "checks passed: $nok, failed: $nbad\n";
print "  $cat{$_}  $_\n" for sort keys %cat;
