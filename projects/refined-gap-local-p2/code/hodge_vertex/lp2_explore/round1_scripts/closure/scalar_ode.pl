use strict; use warnings; use Math::BigRat;
# closure.pl -- closure-lemma test (pathway stage 2, test 1). Lines 1-158 are copied from hodge_vertex/lp2_explore/edges_sing.pl
# (series machinery, rows P0,P1,P2 at c1 = 1 and at the section c1), with ONE change: off shell the section c1 is chosen CONSISTENTLY,
# theta log c1 = XX (shifted), i.e. c1 = C1cf * exp(Int 3 dXX/(s-1) ds), so that the rows are an honest gauge of a section c1 that
# fails only the Picard-Fuchs/Riccati equation. Usage: perl closure.pl cPP_z9.m e1 e2 [NS]; env XXSHIFT=q, XXSHIFT1=q, XXZERO=1 (XX:=0, c1:=1).
# edges_sing.pl -- CHECK 2 (analogue of F3/L2/L3 of hodge_vertex_conifold/leaf_proof/leaf_structure_proof.md) for K_P2 x C^2 on T_F,
# frame B (rows times h_w(z)^-1, w = (E, -e1, -e2), as hv_sample.m), mode-B Hodge vertices, conifold polarisation.
# Exact arithmetic (Math::BigRat, and Q(omega), omega = exp(2 pi i/3), for the two regular colours). No Mathematica.
#
# Construction (from lp2_lib.m buildRmatrix/evalPoint, re-derived by hand; see notes below):
#  * cPP_k (cPP_z9.m) = 3^(p_k/2) sq^-1 N_k(L)/(E+3L)^(k+1), deg N_k <= k (N_k found by exact interpolation, checked).
#  * Normalised rows Phat_{a,i} := sqD_i P_{a,i} (sqD = sqrtDelta of lp2_lib.m), so that the per-half-edge vertex factors are absorbed:
#    Ehat_ij := sqD_i sqD_j Ed_ij, Dilhat_i(n) := sqD_i Dil_i(n), and a coloured graph = prod_v Delta^(g_v-1) * Hodge * prod Ehat * prod Dilhat,
#    exactly the conifold convention of leaf_structure_proof.md s.1.2 (C6, C10, C11).
#    Phat0_k = 3^((3+p_k)/2-k) alpha^2 N_k(L)/(alpha+L)^k,  L = -alpha/(1+3t),  alpha = E/3.
#    Phat1 = L Phat0 + z Dop(Phat0),  Phat2 = Del [L Phat1 + z Dop(Phat1) - z XX Phat1]      (at c1 = 1, i.e. C1i = 1, C1p = XX/y)
#    Dop(F) = (t/3) dF/dt - F/3 + (t/(1+3t)) F,   Del = 1 + 27 y.
#    The edge factors depend on (C1i, C1p) only through XX (the C1i^2 terms cancel: etainv11 = E etainv12), so c1 = 1 is no loss for edges.
#    The true rows for a section c1 are Phat1_true = Phat1/c1, Phat2_true = c1 Phat2 + (E/2)(c1 - 1/c1) Phat1 (reported for c1 = -2F1(1/3,2/3;1;Del),
#    the normalisation C1cf of lp2_lib.m tcSetup: C1cf = 27 y dtc/dDel = -wc).
#  * Colour i at t_i = omega^i t, t real; expansion variable s := 1 + 3t (real branch), so Del = 3s - 3s^2 + s^3 and t_c = Del + O(Del^2) = 3s + O(s^2).
#    (t/3) d/dt at every colour = ((s-1)/3) d/ds.  XX = (Del - 1) wc'(Del)/wc(Del), wc = 2F1(1/3,2/3;1;Del) (the analytic solution of
#    tcSetup's L2b at Del = 0; exponents 0,0, so it is unique; XX(0) = -2/9).
#  * Edges: pre[p][q] = -(-1)^(p+q) Sum_ab etainv_ab Phat_a,i[p] Phat_b,j[q], Ehat_ij(m,n) = Sum_{k=0}^n (-1)^k pre[m+1+k][n-k] (evalPoint).
# Leaf tests (variable psi = z of the edge generating function Sum E(m,n) z^m w^n, as conifold C12; leaves lam_j = 1/(E - j z)):
#   L3-analogue: [s^k] Ehat_00(z,w) in span{lam_i(z) lam_l(w): (i-1)+(l-1) <= k};  L2-analogue: [s^k] [w^b] Ehat_0j(z,w) in span{lam_1..lam_(k+1)}(z).
#   Because s = t_c * unit, and the statements are filtrations by order, they are equivalent for s and t_c, and for rho^e times the factors.
# Usage: perl edges_sing.pl cPP_z9.m e1 e2 [NS=6] [KS=3]      (NS = s-precision of the rows; edges are then known to s^(NS-3))
#   Runs of 2026-11-01: NS = 7, KS = 4 at (e1,e2) = (7,3), (5,2), (4,9) -> edges_sing_out_<e1>_<e2>.txt (about 1 min each).
#   Diagnostics (off-shell polarisation, not part of the check): XXSHIFT=q adds the constant q to XX, XXSHIFT1=q adds q*s.
#   Results: q = 1/10 constant shift -> E_00 acquires an s^-1 pole; q = 1/10 linear shift -> edges stay regular but the rank-one
#   limit and the leaf supports fail. So the checks below are sensitive to the polarisation and pass only on shell.
# Expected output: symplecticity exact (all zeros), all nine edge factors regular at s = 0, frame-B identity row = E/(E+z) at s = 0,
#   E_00 at s^0 of rank one in lam_1 lam_1, E_0j at s^0 = lam_1(z) chi_j(w), leaf supports as in the conifold (if LC_P2 holds).

my ($file, $E1, $E2, $NS, $KS) = @ARGV; $NS //= 6; $KS //= 3;
sub QQ { Math::BigRat->new("$_[0]") }
my $ZERO = QQ(0); my $ONE = QQ(1);
local $/; open my $fh,'<',$file or die "$file: $!"; my $s=<$fh>; close $fh;
$s =~ s/\s+//g; $s =~ s/^\{//; $s =~ s/\}$//;
my @parts; { my $d=0; my $cur='';
for my $ch (split //,$s){ if($ch eq '['||$ch eq '('){$d++} elsif($ch eq ']'||$ch eq ')'){$d--}
  if($ch eq ',' && $d==0){push @parts,$cur;$cur=''} else {$cur.=$ch} } push @parts,$cur; }
my $NZ = $#parts;   # z^0..z^NZ
my $bs = chr(92);
our ($e1,$e2,$s3,$sqv,$L);
$e1 = QQ($E1); $e2 = QQ($E2); my $E = $e1+$e2; my $alpha = $E/3;
sub ev { my ($x,$Lv,$s3v,$sqq)=@_; my $t=$x; $t =~ s/\Q$bs\E//g;
  $t =~ s/\Q[Epsilon][1]\E/\$e1/g; $t =~ s/\Q[Epsilon][2]\E/\$e2/g; $t =~ s/\QLLL[0][y]\E/\$L/g;
  $t =~ s/\QSqrt[3]\E/\$s3/g; $t =~ s/\bsq\b/\$sqv/g; $t =~ s/\^/**/g; $t =~ tr/[]/()/;
  $t =~ s/(?<![\$\w.])(\d+)/QQ('$1')/g;
  local $L = $Lv; local $s3 = $s3v; local $sqv = $sqq;
  my $v = eval $t; die "$@ in $t" if $@; return $v; }
print "e1 = $e1, e2 = $e2, E = $E; z-orders 0..$NZ; s-precision $NS; leaf tests k = 0..$KS\n";
# ---- N_k, p_k ----
my (@pk, @Nc);
for my $k (0..$NZ) {
  my $a = ev($parts[$k], QQ(5), $ONE, $ONE); my $b = ev($parts[$k], QQ(5), QQ(3), $ONE); my $c = ev($parts[$k], QQ(5), $ONE, QQ(2));
  die "sq power not -1 at k=$k" unless $c*2 == $a;
  my $q = $b/$a; my $p = ($q == 3) ? 1 : ($q == QQ('1/3') ? -1 : die "sqrt3 power ? at k=$k");
  push @pk, $p;
  my @xs = map { QQ($_+2) } 0..($k+3); my @ys = map { ev($parts[$k], $_, $ONE, $ONE) * ($E + 3*$_)**($k+1) } @xs;
  my @cc = @ys; my $n = scalar(@xs);
  for my $j (1..$n-1) { for (my $i=$n-1; $i>=$j; $i--) { $cc[$i] = ($cc[$i]-$cc[$i-1])/($xs[$i]-$xs[$i-$j]); } }
  my @poly = ($ZERO) x $n; my @basis = ($ONE);
  for my $i (0..$n-1) { for my $m (0..$#basis) { $poly[$m] += $cc[$i]*$basis[$m]; }
    my @nb = ($ZERO) x (scalar(@basis)+1); for my $m (0..$#basis) { $nb[$m+1] += $basis[$m]; $nb[$m] -= $xs[$i]*$basis[$m]; } @basis = @nb; }
  for my $m ($k+1..$n-1) { die "deg N_$k > $k" if $poly[$m] != 0; }
  $Nc[$k] = [ @poly[0..$k] ];
}
print "Sqrt[3] powers p_k: @pk ; deg N_k <= k and sq^-1 checked for all k\n";
# ---- Q(omega) ----
sub K { [ $_[0], $_[1] // $ZERO ] }
sub kq { [ QQ($_[0]), $ZERO ] }
sub kadd { [ $_[0][0]+$_[1][0], $_[0][1]+$_[1][1] ] }
sub ksub { [ $_[0][0]-$_[1][0], $_[0][1]-$_[1][1] ] }
sub kmul { my ($x,$y)=@_; my $bd = $x->[1]*$y->[1]; [ $x->[0]*$y->[0]-$bd, $x->[0]*$y->[1]+$x->[1]*$y->[0]-$bd ] }
sub ksc { [ $_[0]*$_[1][0], $_[0]*$_[1][1] ] }
sub kz { $_[0][0] == 0 && $_[0][1] == 0 }
sub kinv { my ($x)=@_; my ($a,$b)=@$x; my $nn = $a*$a-$a*$b+$b*$b; die "kinv 0" if $nn == 0; [ ($a-$b)/$nn, -$b/$nn ] }
sub kstr { my ($x)=@_; return "$x->[0]" if $x->[1] == 0; return "($x->[0] + ($x->[1]) w)"; }
my $W = [ $ZERO, $ONE ]; my $W2 = [ -$ONE, -$ONE ];
# ---- truncated Laurent series in s: {v, N, a} (a->[i] = coefficient of s^(v+i)); known for exponents < N; all N capped at NS ----
my $K0 = K($ZERO);
sub sget { my ($f,$e)=@_; return $K0 if $e < $f->{v}; die "precision: s^$e of series known to $f->{N}" if $e >= $f->{N}; return $f->{a}[$e-$f->{v}] // $K0; }
sub snew { my ($v,$N,$a)=@_; $N = $NS if $N > $NS; { v=>$v, N=>$N, a=>$a } }
my $INF = 1000000;
sub sconst { { v=>0, N=>$INF, a=>[ $_[0] ] } }
sub spoly { my ($v,$a)=@_; { v=>$v, N=>$INF, a=>[ @$a ] } }   # exact (Laurent) polynomial s^v * (a0 + a1 s + ...)
sub sext { my ($f)=@_; my $x = $f->{v} + scalar(@{$f->{a}}); return $x < $f->{N} ? $x : $f->{N}; }   # exponents >= sext are zero or unknown
sub sadd { my ($f,$g)=@_; my $v = $f->{v} < $g->{v} ? $f->{v} : $g->{v}; my $N = $f->{N} < $g->{N} ? $f->{N} : $g->{N};
  my $top = $N; if ($N >= $INF) { my ($x,$y) = (sext($f), sext($g)); $top = $x > $y ? $x : $y; }
  my $r = { v=>$v, N=>$N, a=>[ map { kadd(sget($f,$_),sget($g,$_)) } $v..$top-1 ] }; $r->{N} = $NS if $r->{N} > $NS && $r->{N} < $INF; $r }
sub ssub { my ($f,$g)=@_; sadd($f, sscal(K(-$ONE),$g)) }
sub sscal { my ($k,$f)=@_; return { v=>$f->{v}, N=>$f->{N}, a=>[ map { kmul($k,$_ // $K0) } @{$f->{a}} ] }; snew($f->{v},$f->{N},[ map { kmul($k,$_ // $K0) } @{$f->{a}} ]) }
sub smul { my ($f,$g)=@_; my $v = $f->{v}+$g->{v}; my $N1 = $f->{N}+$g->{v}; my $N2 = $g->{N}+$f->{v}; my $N = $N1<$N2?$N1:$N2; $N = $NS if $N > $NS;
  my @a; for my $e ($v..$N-1) { my $acc = $K0; for my $i ($f->{v}..$e-$g->{v}) { my $x = $f->{a}[$i-$f->{v}]; next unless defined $x; next if kz($x);
      my $y = $g->{a}[$e-$i-$g->{v}]; next unless defined $y; $acc = kadd($acc, kmul($x,$y)); } push @a,$acc; }
  { v=>$v, N=>$N, a=>\@a } }
sub sstrip { my ($f)=@_; my %g = (v=>$f->{v}, N=>$f->{N}, a=>[ @{$f->{a}} ]); while ($g{v} < $g{N} && kz($g{a}[0] // $K0)) { shift @{$g{a}}; $g{v}++; } \%g }
sub sinv { my ($f)=@_; my $g = sstrip($f); die "sinv: zero series" if $g->{v} >= $g->{N}; my $R = $g->{N}-$g->{v}; my $Rc = $NS + $g->{v}; $R = $Rc if $R > $Rc;
  my $b0 = kinv($g->{a}[0]); my @b = ($b0);
  for my $n (1..$R-1) { my $acc = $K0; for my $j (1..$n) { my $x = $g->{a}[$j] // $K0; $acc = kadd($acc, kmul($x,$b[$n-$j])); } push @b, kmul(K(-$ONE), kmul($b0,$acc)); }
  snew(-$g->{v}, -$g->{v}+$R, \@b) }
sub sder { my ($f)=@_; my @a; for my $e ($f->{v}..sext($f)-1) { push @a, ksc(QQ($e), sget($f,$e)); }
  return { v=>$f->{v}-1, N=>($f->{N} >= $INF ? $INF : $f->{N}-1), a=>[ @a ] } }   # coefficient e*a_e at s^(e-1)
sub sval { my ($f)=@_; my $g = sstrip($f); return $g->{v} < $g->{N} ? $g->{v} : "none(<$g->{N})"; }
# ---- XX on the conifold section, as a series in s ----
my $sS = spoly(1, [ K($ONE) ]);
my $Del = spoly(1, [ K(QQ(3)), K(QQ(-3)), K($ONE) ]);
my @wa = ($ONE); for my $n (1..$NS+2) { push @wa, $wa[-1] * (QQ($n-1)+QQ('1/3'))*(QQ($n-1)+QQ('2/3'))/(QQ($n)*QQ($n)); }
my $wc = sconst(K($ONE)); my $wcp = sconst(K($wa[1])); { my $Dp = sconst(K($ONE));
  for my $n (1..$NS+1) { $Dp = smul($Dp,$Del); $wc = sadd($wc, sscal(K($wa[$n]),$Dp)); $wcp = sadd($wcp, sscal(K($wa[$n+1]*($n+1)),$Dp)); } }
my $XX = smul(smul(sadd($Del, sconst(K(-$ONE))), $wcp), sinv($wc));
if ($ENV{XXSHIFT}) { $XX = sadd($XX, sconst(K(QQ($ENV{XXSHIFT})))); print "OFF-SHELL TEST: XX shifted by $ENV{XXSHIFT}
"; }   # diagnostic only
if ($ENV{XXSHIFT1}) { $XX = sadd($XX, spoly(1, [ K(QQ($ENV{XXSHIFT1})) ])); print "OFF-SHELL TEST: XX shifted by ($ENV{XXSHIFT1}) s
"; }   # diagnostic only
print "XX(s) = ", join(', ', map { kstr(sget($XX,$_)) } 0..$NS-1), "   (XX(0) must be -2/9)\n";
my $XXon = $XX;   # (possibly shifted) XX actually used in P2
my $dXX = ssub($XX, smul(smul(sadd($Del, sconst(K(-$ONE))), $wcp), sinv($wc)));   # shift
my $c1sec;
if ($ENV{XXZERO}) { $XX = sconst($K0); $XXon = $XX; $c1sec = sconst(K($ONE)); print "OFF-SHELL TEST: XX := 0, c1 := 1 (constant section)\n"; }
else {
  # g = Int_0^s 3 dXX/(s-1) ds ; c1 = -wc * exp(g)
  my $integrand = sscal(K(QQ(3)), smul($dXX, sinv(sadd($sS, sconst(K(-$ONE))))));
  my @ga = ($K0); for my $e (0..$NS-1) { push @ga, ksc($ONE/QQ($e+1), sget($integrand,$e)); }
  my $g = { v=>0, N=>$NS, a=>\@ga };
  my $ex = sconst(K($ONE)); my $gp = sconst(K($ONE)); my $fact = $ONE;
  for my $n (1..$NS) { $gp = smul($gp,$g); $fact *= $n; $ex = sadd($ex, sscal(K($ONE/$fact), $gp)); }
  $c1sec = sscal(K(-$ONE), smul($wc, $ex));
}
# check theta log c1 = XX for the chosen section
my $thetaS = sub { my ($F)=@_; smul(sscal(K(QQ('1/3')), sadd($sS, sconst(K(-$ONE)))), sder($F)) };
{ my $d = sstrip(ssub(smul($thetaS->($c1sec), sinv($c1sec)), $XX));
  print "section check theta log c1 - XX: valuation ", ($d->{v} < $d->{N} ? "$d->{v} (NONZERO)" : "none (zero to s^".($d->{N}-1).")"), "\n"; }
my (@DopC, @LC);
# ---- h_w^-1 (frame B) ----
my @B = (QQ('1/6'),QQ('-1/30'),QQ('1/42'),QQ('-1/30'),QQ('5/66'));
my @w = ($E, -$e1, -$e2); my @lg = ($ZERO) x ($NZ+1);
for my $k (1..5) { my $m = 2*$k-1; next if $m > $NZ; my $S = $ZERO; $S += $ONE/($_**$m) for @w; $lg[$m] = $B[$k-1]/(2*$k*(2*$k-1))*$S; }
my @hinv = ($ONE, ($ZERO) x $NZ);
for my $n (1..$NZ) { my $v = $ZERO; for my $k (1..$n) { $v += $k*$lg[$k]*$hinv[$n-$k]; } $hinv[$n] = $v/$n; }
# ---- rows per colour ----
my @om = (K($ONE), $W, $W2);
my (@rows, @Delta, @rowsTrue);
for my $i (0..2) {
  my $oi = $om[$i]; my $ti = spoly(0, [ ksc(QQ('-1/3'),$oi), ksc(QQ('1/3'),$oi) ]);        # t_i = omega^i (s-1)/3 (exact)
  my ($vi, $ivi);
  if ($i == 0) { $vi = spoly(1, [ K($ONE) ]); $ivi = spoly(-1, [ K($ONE) ]); }            # v_0 = 1+3t_0 = s exactly
  else { $vi = spoly(0, [ ksub(K($ONE),$oi), $oi ]); $ivi = sinv($vi); }
  my $Li = sscal(K(-$alpha), $ivi);
  my $it = sinv($ti);
  # (alpha+L)^-1 = v/(3 alpha t): exact in terms of v and 1/t (no Laurent cancellation needed)
  my $iApL = sscal(K($ONE/(3*$alpha)), smul($vi, $it));
  my $tov = smul($ti, $ivi);
  my $Dop = sub { my ($F)=@_; sadd(sadd(smul(sscal(K(QQ('1/3'))
            , sadd($sS, sconst(K(-$ONE)))), sder($F)), sscal(K(QQ('-1/3')),$F)), smul($tov,$F)) };
  # L^j (alpha+L)^-k = (-alpha)^j v^(k-j) (3 alpha t)^-k : every term regular, no precision loss
  my @Vp = (sconst(K($ONE))); push @Vp, smul($Vp[-1],$vi) for 1..$NZ;
  my @Tp = (sconst(K($ONE))); push @Tp, smul($Tp[-1],$it) for 1..$NZ;
  my @P0; for my $k (0..$NZ) { my $pp = (3+$pk[$k])/2 - $k; my $pref = $pp >= 0 ? QQ(3)**$pp : $ONE/(QQ(3)**(-$pp)); $pref *= $alpha*$alpha/((3*$alpha)**$k);
    my $acc = sconst($K0); for my $j (0..$k) { $acc = sadd($acc, sscal(K($Nc[$k][$j]*(-$alpha)**$j), $Vp[$k-$j])); }
    push @P0, sscal(K($pref), smul($acc, $Tp[$k])); }
  my @P1; for my $m (0..$NZ) { my $x = smul($Li,$P0[$m]); $x = sadd($x, $Dop->($P0[$m-1])) if $m >= 1; push @P1, $x; }
  my @P2; for my $m (0..$NZ) { my $x = smul($Li,$P1[$m]); $x = sadd($x, ssub($Dop->($P1[$m-1]), smul($XX,$P1[$m-1]))) if $m >= 1; push @P2, smul($Del,$x); }
  # frame B
  my @R; for my $a ([@P0],[@P1],[@P2]) { my @r; for my $k (0..$NZ) { my $acc = sconst($K0); for my $j (0..$k) { $acc = sadd($acc, sscal(K($hinv[$k-$j]), $a->[$j])); } push @r, $acc; } push @R, \@r; }
  $rows[$i] = \@R; $DopC[$i] = $Dop; $LC[$i] = $Li;
  $Delta[$i] = sscal(K(3*$e1*$e2*$E**3), smul(smul($ti,$ti), smul($ivi,$ivi)));
  my $ic1 = sinv($c1sec);
  $rowsTrue[$i] = [ $R[0], [ map { smul($R[1][$_], $ic1) } 0..$NZ ],
     [ map { sadd(smul($c1sec,$R[2][$_]), smul(sscal(K($E/2), ssub($c1sec,$ic1)), $R[1][$_])) } 0..$NZ ] ];
}
print "rows built\n";
# ---- scalar check: equivariant Picard-Fuchs (BS eq:pf on T_F) in the e^{u/z}-stripped form ----
# Del Th^3 P0 + (E + 27 y z) Th^2 P0 + (E^2/3 + 6 y z^2) Th P0 + E^3/27 P0 = 0,  Th = L + z Dop (sqD-conjugated, frame B rows)
my $y27 = sadd($Del, sconst(K(-$ONE)));   # 27 y
for my $i (0..2) { my $P = $rows[$i][0];
  my $Th = sub { my ($F)=@_; [ map { my $m = $_; my $x = smul($LC[$i], $F->[$m]); $x = sadd($x, $DopC[$i]->($F->[$m-1])) if $m >= 1; $x } 0..$NZ ] };
  my $T1 = $Th->($P); my $T2 = $Th->($T1); my $T3 = $Th->($T2);
  my $ok = 1; my $mn = 99;
  for my $m (0..$NZ) { my $acc = smul($Del, $T3->[$m]);
    $acc = sadd($acc, sscal(K($E), $T2->[$m])); $acc = sadd($acc, smul($y27, $T2->[$m-1])) if $m >= 1;
    $acc = sadd($acc, sscal(K($E*$E/3), $T1->[$m])); $acc = sadd($acc, smul(sscal(K(QQ($ENV{PFCTRL} ? '1/9' : '2/9')),$y27), $T1->[$m-2])) if $m >= 2;
    $acc = sadd($acc, sscal(K($E**3/27), $P->[$m]));
    my $g = sstrip($acc); $ok = 0 if $g->{v} < $g->{N}; $mn = $g->{N} if $g->{N} < $mn; }
  print "colour $i: PF operator annihilates Phat0 (frame B), z^0..z^$NZ: ", ($ok ? "OK" : "FAILED"), " (exact to s^", $mn-1, ")\n"; }
print "done\n";
