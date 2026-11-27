# h2.pl -- Lemma H2 test (pathway stage 3, test 2). Lines 1-154 of lp2_explore/edges_sing.pl copied verbatim (row construction),
# restricted to colour 0; the H2 analysis is appended below. Usage: perl h2.pl cPP_z9.m e1 e2 NS [KMAX]   (env XXSHIFT, XXSHIFT1 as edges_sing.pl)
use strict; use warnings; use Math::BigRat;
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
my $c1sec = sscal(K(-$ONE), $wc);   # C1cf = -wc
# ---- h_w^-1 (frame B) ----
my @B = (QQ('1/6'),QQ('-1/30'),QQ('1/42'),QQ('-1/30'),QQ('5/66'));
my @w = ($E, -$e1, -$e2); my @lg = ($ZERO) x ($NZ+1);
for my $k (1..5) { my $m = 2*$k-1; next if $m > $NZ; my $S = $ZERO; $S += $ONE/($_**$m) for @w; $lg[$m] = $B[$k-1]/(2*$k*(2*$k-1))*$S; }
my @hinv = ($ONE, ($ZERO) x $NZ);
for my $n (1..$NZ) { my $v = $ZERO; for my $k (1..$n) { $v += $k*$lg[$k]*$hinv[$n-$k]; } $hinv[$n] = $v/$n; }
# ---- rows per colour ----
my @om = (K($ONE), $W, $W2);
my (@rows, @Delta, @rowsTrue);
for my $i (0..0) {
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
  $rows[$i] = \@R;
  $Delta[$i] = sscal(K(3*$e1*$e2*$E**3), smul(smul($ti,$ti), smul($ivi,$ivi)));
  my $ic1 = sinv($c1sec);
  $rowsTrue[$i] = [ $R[0], [ map { smul($R[1][$_], $ic1) } 0..$NZ ],
     [ map { sadd(smul($c1sec,$R[2][$_]), smul(sscal(K($E/2), ssub($c1sec,$ic1)), $R[1][$_])) } 0..$NZ ] ];
}
print "rows built\n";
# ===================== H2 analysis (new code) =====================
# sigma(psi) := rho * (Phat0, Phat1_true, Phat2_true)(z = -psi), frame B, colour 0, with rho := -s/(t t_c) (rho -> 1 at s = 0;
# t_c^2 rho^2 Delta_0 = 3 e1 e2 E^3 exactly), true section c1 = C1cf = -wc. Also the c1 = 1 code rows (gauge check).
# For each t_c-order k: F(psi) = [t_c^k] sigma_a; find minimal J with F prod_{j<=J}(E - j psi) a polynomial of degree <= J+1
# (psi-linear part allowed), then F = a1 psi + a0 + Sum_j beta_j lam_j(psi), lam_j = 1/(E - j psi).
my $KMAX = $KS;
sub re { my ($x)=@_; die "non-real entry at colour 0" unless $x->[1] == 0; $x->[0] }
my $M2 = $NS+3;
my @wcD = (QQ(1)); for my $n (0..$M2) { $wcD[$n+1] = $wcD[$n]*(QQ($n)+QQ('1/3'))*(QQ($n)+QQ('2/3'))/(QQ($n+1)**2); }
my @gg = map { my $n=$_; my $v=QQ(0); $v += $wcD[$_] for 0..$n; $v } 0..$M2;
my @tcD = (QQ(0), map { $gg[$_-1]/$_ } 1..$M2);
my $tc = sconst($K0); { my $Dp = sconst(K($ONE)); for my $n (1..$NS+1) { $Dp = smul($Dp,$Del); $tc = sadd($tc, sscal(K($tcD[$n]),$Dp)); } }
my $t0 = spoly(0,[K(QQ('-1/3')),K(QQ('1/3'))]);
my $rho = sscal(K(-$ONE), smul(smul($sS, sinv($t0)), sinv($tc)));
print "t_c(s) = ", join(', ', map { kstr(sget($tc,$_)) } 1..4), " ...;  rho(s) = ", join(', ', map { kstr(sget($rho,$_)) } 0..3), " ... (known to s^", $rho->{N}-1, ")\n";
{ my $Dl = sscal(K(3*$e1*$e2*$E**3), smul(smul($t0,$t0), smul(spoly(-1,[K($ONE)]),spoly(-1,[K($ONE)]))));
  my $x = sstrip(ssub(smul(smul($tc,$tc),smul(smul($rho,$rho),$Dl)), sconst(K(3*$e1*$e2*$E**3))));
  print "t_c^2 rho^2 Delta_0 - 3 e1 e2 E^3 = 0 to s^", $x->{N}-1, ": ", ($x->{v} >= $x->{N} ? "OK" : "NO"), "\n"; }
my $KT = $KMAX + 3; $KT = $tc->{N}-2 if $KT > $tc->{N}-2;
sub amul { my ($a,$b,$K)=@_; my @c=(QQ(0))x($K+1); for my $i (0..$K){ my $x=$a->[$i]; next unless defined $x; next if $x==0; for my $j (0..$K-$i){ my $y=$b->[$j]; next unless defined $y; $c[$i+$j]+=$x*$y; } } \@c }
sub ainv { my ($a,$K)=@_; my @b=(QQ(1)/$a->[0]); for my $n (1..$K){ my $v=QQ(0); for my $j (1..$n){ $v += ($a->[$j]//QQ(0))*$b[$n-$j]; } $b[$n] = -$v*$b[0]; } \@b }
my @T = map { re(sget($tc,$_+1)) } 0..$KT;
my $sig = [ QQ(1)/$T[0] ];
for my $it (0..$KT+1) { my $x = [ QQ(0), @$sig ]; $#$x = $KT if $#$x > $KT;
  my $comp = [ (QQ(0)) x ($KT+1) ]; my $p = [ QQ(1) ]; for my $n (0..$KT) { for my $i (0..$KT) { $comp->[$i] += $T[$n]*($p->[$i]//QQ(0)); } $p = amul($p,$x,$KT); }
  $sig = ainv($comp,$KT); }
{ my $x = [ QQ(0), @$sig ]; $#$x = $KT; my $comp = [ (QQ(0)) x ($KT+1) ]; my $p = [ QQ(1) ];
  for my $n (0..$KT) { for my $i (0..$KT) { $comp->[$i] += $T[$n]*($p->[$i]//QQ(0)); } $p = amul($p,$x,$KT); }
  my $chk = amul($comp,$sig,$KT); my $ok = ($chk->[0] == 1); $ok &&= ($chk->[$_] == 0) for 1..$KT; print "inverse series s(t_c) check: ", ($ok?"OK":"NO"), "\n"; }
my %sigp; $sigp{0} = [ QQ(1), (QQ(0)) x $KT ]; $sigp{-1} = ainv($sig,$KT); $sigp{-2} = amul($sigp{-1},$sigp{-1},$KT);
for my $e (1..$KT) { $sigp{$e} = amul($sigp{$e-1},$sig,$KT); }
print "s(t_c)/t_c = ", join(', ', map { $sig->[$_] } 0..3), " ...\n";
sub totc { my ($f,$k)=@_;
  die "precision: t_c^$k needs s^$k, known to s^".($f->{N}-1) if $k >= $f->{N};
  die "valuation < -2" if $f->{v} < -2 && grep { !kz($_) } @{$f->{a}}[0..(-3-$f->{v})];
  my $v = QQ(0); for my $e ($f->{v}..$k) { next if $e < -2; my $c = re(sget($f,$e)); next if $c == 0; $v += $c * $sigp{$e}[$k-$e]; } $v }
my $isS = spoly(-1,[K($ONE)]);
my $DopS = sub { my ($F)=@_; sadd(sadd(smul(sscal(K(QQ('1/3')), sadd($sS, sconst(K(-$ONE)))), sder($F)), sscal(K(QQ('-1/3')),$F)), smul(smul($t0,$isS),$F)) };
my $Hinf = sscal(K($E), smul($wc, smul($t0,$isS)));                  # predicted z -> oo limit of the c1=1 row Phat1: E wc t/s
my $Xc = ssub($DopS->($Hinf), smul($XX,$Hinf));                      # Dop(H) - XX H
my $predT = sscal(K(-$ONE), smul(smul($rho, smul($c1sec,$Del)), $Xc)); # psi-coefficient of sigma_2 (true gauge)
my $pred1 = sscal(K(-$ONE), smul(smul($rho, $Del), $Xc));             # psi-coefficient of rho Phat2 (c1 = 1 gauge)
{ my $th = sscal(K(QQ('1/3')), smul(sadd($sS, sconst(K(-$ONE))), smul(sder($wc), sinv($wc))));
  my $x = ssub($XX, $th); my $top = $x->{N}-1; $top = 6 if $top > 6;
  print "identity check on the s-series: XX - (t/3) d/dt log wc = ", join(', ', map { kstr(sget($x,$_)) } 0..$top), " (s^0..s^$top)\n"; }
my %SIG;
for my $g ('true','c1=1') { my $R = $g eq 'true' ? $rowsTrue[0] : $rows[0];
  for my $a (0..2) { for my $m (0..$NZ) { $SIG{$g}[$a][$m] = smul($rho, $R->[$a][$m]); } } }
{ my $minN = 99; for my $g (keys %SIG) { for my $a (0..2) { for my $m (0..$NZ) { $minN = $SIG{$g}[$a][$m]{N} if $SIG{$g}[$a][$m]{N} < $minN; } } }
  print "sigma known to s^", $minN-1, " (so t_c-orders <= ", $minN-1, "); requested KMAX = $KMAX\n"; $KMAX = $minN-1 if $KMAX > $minN-1; }
sub decomp { my ($F,$J,$lin)=@_;
  my @P=(QQ(1)); for my $j (1..$J){ my @n=(QQ(0))x(@P+1); for my $d (0..$#P){ $n[$d]+=$E*$P[$d]; $n[$d+1]-=$j*$P[$d]; } @P=@n; }
  my @N; for my $m (0..$NZ){ my $v=QQ(0); for my $d (0..($m<$#P?$m:$#P)){ $v+=$P[$d]*$F->[$m-$d]; } $N[$m]=$v; }
  my $deg = $J + ($lin?1:0); return { untestable=>1 } if $deg > $NZ;
  for my $i ($deg+1..$NZ) { return undef if $N[$i] != 0; }
  my @r = @N[0..$deg]; my @q = (QQ(0)) x ($deg-$J+1);
  for (my $d=$deg; $d>=$J; $d--) { my $c = $r[$d]/$P[$J]; $q[$d-$J] = $c; for my $i (0..$J) { $r[$d-$J+$i] -= $c*$P[$i]; } }
  my @beta; for my $j (1..$J) { my $x = $E/$j; my $v=QQ(0); for my $i (reverse 0..$J-1) { $v = $v*$x + $r[$i]; }
    my $den = $ONE; for my $i (1..$J) { next if $i==$j; $den *= ($E - $i*$x); } $beta[$j] = $v/$den; }
  return { a0=>$q[0], a1=>($lin ? $q[1] : QQ(0)), beta=>\@beta, J=>$J, checks=>$NZ-$deg }; }
sub fmtd { my ($d)=@_; my @sup = grep { $d->{beta}[$_] != 0 } 1..$d->{J}; my $lim = QQ(0); $lim += $d->{beta}[$_]/$_ for 1..$d->{J};
  return "a1(psi) = $d->{a1}; a0(lam_0) = $d->{a0}; leaves {".join(',',@sup)."}; Sum beta_j/j = $lim"; }
print "\n";
my %bad;
for my $g ('true','c1=1') {
  print "=== gauge $g ===\n";
  for my $k (-1..$KMAX) {
    for my $a (0..2) {
      my @F = map { (($_ % 2) ? -1 : 1) * totc($SIG{$g}[$a][$_], $k) } 0..$NZ;
      my $dm; for my $J (0..$NZ) { my $d = decomp(\@F,$J,1); next unless $d; next if $d->{untestable}; next if $d->{checks} < 1; $dm = $d; last; }
      my $df = decomp(\@F,$k+1 < 0 ? 0 : $k+1,1); my $dfu = ($df && $df->{untestable}) ? 1 : 0; $df = undef if $dfu;
      my $str = sprintf "  t_c^%-2d sigma_%d: ", $k, $a;
      if ($dm) { $str .= "minimal J = $dm->{J} ($dm->{checks} checks): ".fmtd($dm); } else { $str .= "minimal J: untestable at z^$NZ"; }
      if ($df && (!$dm || $dm->{J} != $df->{J})) { $str .= "\n      [assuming pole set j <= ".($k+1).": ".($df->{checks})." checks] ".fmtd($df); }
      if (!$df && !$dfu) { $str .= "\n      [pole set j <= ".($k+1)." with psi-linear term REJECTED by the data]"; }
      if ($dfu && !$dm) { $str .= "\n      [pole set j <= ".($k+1)." plus psi-linear term: too few z-orders]"; }
      if (!$dm) { my $d0; for my $J (0..$NZ) { my $d = decomp(\@F,$J,0); next unless $d; next if $d->{untestable}; next if $d->{checks} < 1; $d0 = $d; last; }
        $str .= "\n      [hypothesis NO psi-linear term: ".($d0 ? "minimal J = $d0->{J} ($d0->{checks} checks): ".fmtd($d0) : "untestable")."]"; }
      if ($a == 2) { my $p = totc($g eq 'true' ? $predT : $pred1, $k); $str .= "\n      closed-form predicted psi-coefficient: $p";
        my $d = $dm // $df; $str .= ($d ? ($d->{a1} == $p ? "  (= fitted)" : "  (DIFFERS from fitted $d->{a1})") : ""); }
      print "$str\n";
      my $d = $dm // $df;
      if ($k >= 0) { $bad{"$g k=$k a=$a: J>k+1"}=1 if ($dm && $dm->{J} > $k+1) || (!$dm && !$df && !$dfu); $bad{"$g k=$k a=$a: psi-linear"}=1 if $d && $d->{a1} != 0; }
      if ($k == -1 && $d) { $bad{"$g sigma^(-1)_$a psi-dependent"}=1 if $d->{a1} != 0 || grep { $d->{beta}[$_] != 0 } 1..$d->{J}; }
    }
    if ($k >= 0 && $k <= 2 && $g eq 'true') { my @F = map { (($_ % 2) ? -1 : 1) * totc($SIG{$g}[2][$_], $k) } 0..$NZ; my $d = decomp(\@F,$k+1,1);
      if ($d) { print "      sigma_2^($k) beta_j (j=1..",$k+1,"): ", join(', ', map { $d->{beta}[$_] } 1..$k+1), "\n"; } }
  }
}
# ---- structural identities (true gauge), exact as s-series, every z-order:
#   (I1) sigma_1 = Dc G,  (I2) sigma_2 = kappa Dc sigma_1 + mu sigma_1,   G := sigma_0 = rho Phat0 (in z),
#   Dc F := (L_1 + z/t_c) F + z d/dt_c F,  d/dt_c = (1/c1)(t/3) d/dt = (1/c1)((s-1)/3) d/ds,  L_1 := -E/(3 s c1),
#   kappa := c1^3 Del,  mu := (E/2)(c1^2 - 1).  (I2) uses XX = c1 d/dt_c log c1 (= (t/3) d/dt log wc); off shell it fails.
{ my $ic1 = sinv($c1sec); my $itc = sinv($tc);
  my $dc = sub { smul($ic1, smul(sscal(K(QQ('1/3')), sadd($sS, sconst(K(-$ONE)))), sder($_[0]))) };
  my $L1 = sscal(K(-$E/3), smul($isS, $ic1));
  my $kap = smul(smul($c1sec, smul($c1sec,$c1sec)), $Del); my $mu = sscal(K($E/2), sadd(smul($c1sec,$c1sec), sconst(K(-$ONE))));
  my $Dc = sub { my ($F,$m)=@_; my $x = smul($L1, $F->[$m]); $x = sadd($x, sadd(smul($itc,$F->[$m-1]), $dc->($F->[$m-1]))) if $m >= 1; $x };
  my @G = @{$SIG{true}[0]}; my @S1 = @{$SIG{true}[1]}; my @S2 = @{$SIG{true}[2]};
  my ($ok1,$ok2,$n1,$n2) = (1,1,99,99); my @res2;
  for my $m (0..$NZ) { my $d1 = sstrip(ssub($S1[$m], $Dc->(\@G,$m))); $n1 = $d1->{N} if $d1->{N} < $n1; $ok1 = 0 if $d1->{v} < $d1->{N};
    my $d2 = sstrip(ssub($S2[$m], sadd(smul($kap, $Dc->(\@S1,$m)), smul($mu,$S1[$m])))); $n2 = $d2->{N} if $d2->{N} < $n2; $ok2 = 0 if $d2->{v} < $d2->{N}; push @res2, sval($d2); }
  print "\n(I1) sigma_1 = Dc sigma_0 exactly (all z^0..z^$NZ, to s^", $n1-1, "): ", ($ok1?"OK":"NO"), "\n";
  print "(I2) sigma_2 = c1^3 Del Dc sigma_1 + (E/2)(c1^2-1) sigma_1 exactly (to s^", $n2-1, "): ", ($ok2?"OK":"NO (residual s-valuations per z-order: @res2)"), "\n";
  print "    L_1 = -E/(3 s c1): t_c-coefficients t_c^-1..t_c^3: ", join(', ', map { totc($L1,$_) } -1..3), "\n";
  print "    kappa = c1^3 Del: t_c^0..t_c^3: ", join(', ', map { totc($kap,$_) } 0..3), ";  mu: t_c^0..t_c^3: ", join(', ', map { totc($mu,$_) } 0..3), "\n";
}
print "\nVIOLATIONS of (sigma^(-1) psi-independent; sigma^(k) in (Lam^+_{k+1})^3, no psi-linear term): ", (%bad ? join('; ', sort keys %bad) : "none"), "\n";
print "done\n";
