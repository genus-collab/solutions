use strict; use warnings; use Math::BigRat;
# hrow.pl -- CHECK 1 of the local-P2 regularity scoping (2026-11-01): leaf support of the frame-B H-row of K_P2 x C^2 (T_F)
# at the singular idempotent near the conifold, exact (Math::BigRat) at rational eps1, eps2.
#
# Definitions (cluster_genus4_hodge/lp2_lib.m, buildRmatrix): P0 = Sum_k cPP_k z^k (identity column, cPP_z9.m, LLL[0][y] -> L),
#   P1 = C1i (L P0 + (z t/3) dP0/dt),  L = -alpha/(1+3t), alpha = E/3, E = eps1+eps2, u := 1/L = -(1+3t)/alpha.
# Facts used/checked: cPP_0 = 1/sqDsym exactly (sqDsym of lp2_lib.m), so r := P0/cPP_0 = sqD P0 is the normalised identity row
#   (r[z^0] = 1) and N1 := P1/(C1i cPP_0) = sqD P1/C1i is the normalised H-row. In u:
#       N1 = ((E+z)/E) r/u + (z/E)(1 + E u/3) dr/du                                         (*)
#   so [u^-1] N1 = (E+z)/E r(z;0) and [u^m] N1 = ((E+(m+2)z)/E) r_{m+1} + (m z/3) r_m     (m >= 0), r_m := [u^m] r.
#   Frame B: every row times h_w(z)^-1, w = (E, -e1, -e2) (t-independent, so (*) holds verbatim in frame B).
# The script computes N1 twice: directly from the u-series of cPP_k (no use of (*)) and from (*), and checks equality.
# Tests, for every u-order m = -1..MU and z^0..z^9 (the z-truncation of cPP_z9.m):
#   (i)   [u^-1] N1B = 1 (z-independent: analogue of the conifold c/t);
#   (ii)  minimal J with [u^m] N1B * prod_{j=1}^J (E + j z) a polynomial of degree <= J (z^(J+1)..z^9 vanish):
#         membership in span{1, 1/(E+jz) : j <= J}; expected J = m+1 (conifold Prop. 6.1), naive bound m+2;
#   (iii) which simple poles z = -E/j occur (residue != 0), and the constant term (value at z = oo);
#   (iv)  the same for the identity row rB = r h^-1 (L1, j <= m+1), as a reference;
#   (v)   the constants with the canonical normalisation rho = 3/(3 + E u) (t^2 rho^2 Delta = const in the u coordinate).
# Usage: perl hrow.pl <cPP_z9.m> [e1 e2 MU]   (default 7 3 6; MU <= 7, since z^9 leaves 9-J vanishing checks).
# Expected: (i) holds; J = m+1 at every m; all poles simple; (vi) in t_c with rho_c and C1i = -1/wc: no constant term at any t_c^m, m >= 0
# (exact analogue of conifold Prop. 6.1, [t^k] sigma_1 in Lam_{k+1}), residue E at t_c^-1 after normalising rho(0) = 1 (conifold: c/t).
# Outputs: hrow_out_<e1>_<e2>.txt.
my $f = shift; my ($E1,$E2) = (shift // 7, shift // 3); my $MU = shift // 6;
local $/; open my $fh,'<',$f or die; my $s=<$fh>; close $fh;
$s =~ s/\s+//g; $s =~ s/^\{//; $s =~ s/\}$//;
my @parts; my $d=0; my $cur='';
for my $ch (split //,$s){ if($ch eq '['||$ch eq '('){$d++} elsif($ch eq ']'||$ch eq ')'){$d--}
  if($ch eq ',' && $d==0){push @parts,$cur;$cur=''} else {$cur.=$ch} } push @parts,$cur;
my $N = $#parts; my $bs = chr(92);
our ($e1,$e2,$s3,$sqv,$L);
sub QQ { Math::BigRat->new(shift) }
$e1 = QQ($E1); $e2 = QQ($E2); my $E = $e1+$e2;
print "eps = ($E1,$E2), E = $E, z-order N = $N, u-orders -1..$MU\n";
sub ev { my ($x,$Lv,$s3v)=@_; my $t=$x; $t =~ s/\Q$bs\E//g;
  $t =~ s/\Q[Epsilon][1]\E/\$e1/g; $t =~ s/\Q[Epsilon][2]\E/\$e2/g; $t =~ s/\QLLL[0][y]\E/\$L/g;
  $t =~ s/\QSqrt[3]\E/\$s3/g; $t =~ s/\bsq\b/\$sqv/g; $t =~ s/\^/**/g; $t =~ tr/[]/()/;
  $t =~ s/(?<![\$\w.])(\d+)/QQ('$1')/g;
  local $L = $Lv; local $s3 = $s3v; local $sqv = QQ(1);
  my $v = eval $t; die "$@ in $t" if $@; return $v; }
my @fac;
for my $k (0..$N) { my $a = ev($parts[$k], QQ(5), QQ(1)); my $b = ev($parts[$k], QQ(5), QQ(3)); my $q = $b/$a;
  push @fac, ($q == 3) ? 1 : ($q == QQ(1)/3 ? -1 : ($q==1?0:die "sqrt3 power?")); }
# cPP_0 * sqD = 1 check (sq -> 1, Sqrt[3] -> s3 symbolic power handled: cPP_0 = Sqrt3/(sq E^2 (E+3L)), sqD = Sqrt3 sq E^2 (E+3L)/3)
{ my $ok = 1; for my $Lv (QQ(2),QQ(5),QQ('7/3')) { my $c0 = ev($parts[0],$Lv,QQ(1)); my $sqD = $E**2*($E+3*$Lv)/3; # times Sqrt3^(+1) sq; c0 carries Sqrt3^(+1)/sq
    $ok = 0 unless $c0*$sqD*3 == 1 && $fac[0] == 1; }   # true value: (c0 Sqrt3)(sqD Sqrt3) = 3 c0 sqD
  print "cPP_0 * sqDsym = 1 exactly at three L values: ", ($ok?"yes":"NO"), "\n"; }
my $MX = $MU + 2;   # u-orders of r needed: 0..MU+1 (and one more for the direct check)
my @Ncoef;
for my $k (0..$N) {
  my @xs = map { QQ($_+2) } 0..($k+2); my @ys = map { ev($parts[$k], $_, QQ(1)) * ($E + 3*$_)**($k+1) } @xs;
  my @c = @ys; my $n = scalar(@xs);
  for my $j (1..$n-1) { for (my $i=$n-1; $i>=$j; $i--) { $c[$i] = ($c[$i]-$c[$i-1])/($xs[$i]-$xs[$i-$j]); } }
  my @poly = (QQ(0)) x $n; my @basis = (QQ(1));
  for my $i (0..$n-1) { for my $m (0..$#basis) { $poly[$m] += $c[$i]*$basis[$m]; }
    my @nb = (QQ(0)) x (scalar(@basis)+1); for my $m (0..$#basis) { $nb[$m+1] += $basis[$m]; $nb[$m] -= $xs[$i]*$basis[$m]; } @basis = @nb; }
  my $deg = 0; for my $m (0..$n-1) { $deg = $m if $poly[$m] != 0; }
  die "k=$k deg N_k = $deg > k" if $deg > $k;
  $Ncoef[$k] = [ @poly[0..$k] ];
}
sub sermul { my ($a,$b,$M)=@_; my @c=(QQ(0))x($M+1); for my $i(0..$M){for my $j(0..$M-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}} return \@c; }
my @inv3Eu = map { (QQ(-1)*$E)**$_ / QQ(3)**($_+1) } 0..$MX+1;   # 1/(3+Eu)
# A_k(u) := cPP_k / (3^(p0/2) sq^-1 n_0) = 3^((p_k-p_0)/2) u [Sum_i n_{k,i} u^(k-i)] / (n_0 (3+Eu)^(k+1));  stored as u-series 0..MX+1
my @A; my @r;
for my $k (0..$N) {
  my $pp = $fac[$k] - $fac[0]; die "odd sqrt3 power" if $pp % 2; my $pref = QQ(3)**($pp/2);
  my @num = (QQ(0))x($MX+2); for my $i (0..$k) { my $e = $k-$i+1; $num[$e] = $Ncoef[$k][$i]*$pref/$Ncoef[0][0] if $e <= $MX+1; }
  my $ser = \@num; for (1..$k+1) { $ser = sermul($ser, \@inv3Eu, $MX+1); }
  $A[$k] = $ser;
  my @n2 = (QQ(0))x($MX+1); for my $i (0..$k) { my $e = $k-$i; $n2[$e] = $Ncoef[$k][$i]*$pref/$Ncoef[0][0] if $e <= $MX; }
  my $s2 = \@n2; for (1..$k) { $s2 = sermul($s2, \@inv3Eu, $MX); } $r[$k] = $s2;   # r_k(u) = cPP_k/cPP_0
}
# h_w^-1
my @B = (QQ('1/6'),QQ('-1/30'),QQ('1/42'),QQ('-1/30'),QQ('5/66'));
my @w = ($E, -$e1, -$e2); my @logh = (QQ(0)) x ($N+1);
for my $k (1..5) { my $m = 2*$k-1; next if $m > $N; my $S = QQ(0); $S += QQ(1)/($_**$m) for @w; $logh[$m] = -$B[$k-1]/(2*$k*(2*$k-1))*$S; }
my @neg = map { -$_ } @logh; my @hinv = (QQ(1), (QQ(0)) x $N);
for my $n (1..$N) { my $v = QQ(0); for my $k (1..$n) { $v += $k*$neg[$k]*$hinv[$n-$k]; } $hinv[$n] = $v/$n; }
sub zmul { my ($a,$b)=@_; my @c=(QQ(0))x($N+1); for my $i(0..$N){for my $j(0..$N-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}} return \@c; }
# z-series of [u^m] X: X given as  coef(k, m)
sub zser { my ($get,$m)=@_; return [ map { $get->($_,$m) } 0..$N ]; }
# ---- direct H-row:  H_k(u) = [z^k] (P1/C1i)/K = A_k/u + ((1+alpha u)/E) dA_{k-1}/du ;  N1 = H / A_0 ----
my $alpha = $E/3;
my @H; for my $k (0..$N) { my @h = (QQ(0))x($MX+1);   # index i <-> u^(i-1), i = 0..MX  (u^-1..u^(MX-1))
  for my $i (0..$MX) { $h[$i] += $A[$k][$i]; }       # A_k/u: [u^(i-1)] = A_k[i]
  if ($k >= 1) { my @dA = map { ($_+1)*$A[$k-1][$_+1] } 0..$MX;   # dA/du, [u^n]
    for my $n (0..$MX) { my $v = $dA[$n]/$E; $v += $alpha*$dA[$n-1]/$E if $n >= 1; $h[$n+1] += $v if $n+1 <= $MX; } }
  $H[$k] = \@h; }
# 1/A_0 as Laurent series: A_0 = u/(3+Eu) * 1 (since num = n_0/n_0 = 1 at e=1)  ->  1/A_0 = (3+Eu)/u
my @N1d; for my $k (0..$N) { my @v = (QQ(0))x($MX);   # index i <-> u^(i-1)
  # (3+Eu)/u * H, H index i <-> u^(i-1)  -> product index: u^(i-2) * 3 + u^(i-1) * E
  for my $i (0..$MX) { my $h = $H[$k][$i]; next if $h == 0;
    my $e3 = $i-2; my $eE = $i-1;   # u-powers
    die "u^-2 term" if $e3 == -2 && $h != 0;
    $v[$e3+1] += 3*$h if $e3+1 >= 0 && $e3+1 < $MX; $v[$eE+1] += $E*$h if $eE+1 < $MX; }
  $N1d[$k] = \@v; }
# formula (*): [u^-1] = (E+z)/E r_0 ; [u^m] = ((E+(m+2)z)/E) r_{m+1} + (m z/3) r_m ; in z-coefficients
sub N1f { my ($k,$m)=@_; # z^k coefficient of [u^m] N1 (point frame)
  if ($m == -1) { my $v = $r[$k][0]; $v += $r[$k-1][0]/$E if $k>=1; return $v; }
  my $v = $r[$k][$m+1]; $v += ($m+2)*$r[$k-1][$m+1]/$E + $m*$r[$k-1][$m]/3 if $k >= 1; return $v; }
my $agree = 1; for my $m (-1..$MU) { for my $k (0..$N) { $agree = 0 unless $N1d[$k][$m+1] == N1f($k,$m); } }
print "direct H-row (from lp2_lib definition) == formula (*) at every u-order -1..$MU, z^0..z^$N: ", ($agree?"yes":"NO"), "\n";
# ---- analysis of a z-series a (z^0..z^N) ----
sub analyse { my ($a,$Jmax)=@_;
  for my $J (0..$Jmax) { my @c = @$a;
    for my $j (1..$J) { my @nc = (QQ(0))x($N+1); for my $i (0..$N) { $nc[$i] += $E*$c[$i]; $nc[$i+1] += $j*$c[$i] if $i < $N; } @c = @nc; }
    my $ok = 1; for my $i ($J+1..$N) { $ok = 0 if $c[$i] != 0; }
    next unless $ok;
    # numerator polynomial c_0..c_J; residues at z = -E/j : nonzero iff numerator(-E/j) != 0
    my @poles; for my $j (1..$J) { my $z0 = -$E/$j; my $v = QQ(0); for my $i (reverse 0..$J) { $v = $v*$z0 + $c[$i]; } push @poles, $j if $v != 0; }
    my $fact = QQ(1); $fact *= $_ for 1..$J; my $const = $c[$J]/$fact;   # value at z = oo
    return ($J, \@poles, $const, $N-$J); }
  return (-1, [], 'n/a', 0); }
sub fmt { my ($lab,$a,$m,$var)=@_; $var //= "u"; my ($J,$p,$cst,$nchk) = analyse($a, $N-1);
  my $exp = $m+1; $exp = 0 if $m == -1;
  printf "  %-22s %s^%-2d: J = %2s (expected <= %d; %d vanishing checks), poles j = {%s}, const(z=oo) = %s\n", $lab, $var, $m, $J, $exp, $nchk, join(',',@$p), $cst;
  return ($J <= $exp && $J >= 0) ? 1 : 0; }
my $allH = 1; my $allR = 1;
print "\nIdentity row rB = r h^-1 (L1 reference; span{1/(E+jz), j <= m+1}):\n";
for my $m (0..$MU) { $allR &= fmt("rB", zmul(zser(sub{ $r[$_[0]][$_[1]] },$m), \@hinv), $m); }
print "\nH-row N1B = N1 h^-1, normalised by cPP_0 = 1/sqD (frame B):\n";
my $res = zmul(zser(\&N1f,-1), \@hinv);
print "  [u^-1] N1B z-coefficients: ", join(', ', @$res), "   (expected 1, 0, ..., 0)\n";
for my $m (0..$MU) { $allH &= fmt("N1B", zmul(zser(\&N1f,$m), \@hinv), $m); }
# canonical rho = 3/(3+Eu) = Sum (-E/3)^n u^n
print "\nH-row with rho = 3/(3+Eu) (rho^2 u^2 Delta = const):\n";
my @rho = map { (QQ(-1)*$E/3)**$_ } 0..$MX;
for my $m (-1..$MU) { my @a = map { my $k=$_; my $v=QQ(0); for my $n (0..$m+1) { $v += $rho[$n]*N1f($k,$m-$n); } $v } 0..$N;
  my $ab = zmul(\@a,\@hinv); if ($m == -1) { print "  [u^-1] rho N1B: ", join(', ', @$ab), "\n"; next; } fmt("rho N1B", $ab, $m); }
print "\nSUMMARY eps=($E1,$E2): identity row J <= m+1 for m = 0..$MU: ", ($allR?"yes":"NO"), "; H-row J <= m+1 for m = 0..$MU: ", ($allH?"yes":"NO"), "\n";
# ---- (vi) conifold flat coordinate t_c and conifold polarisation (hodge_vertex_conifold/lp2_tc_section.m, tcSetup) ----
# Del = 1 + 27 y = 1 + 27 t^3 = -alpha u (3 + 3 alpha u + alpha^2 u^2);  wc(Del) = 2F1(1/3, 2/3; 1; Del) (the solution of L2b with a0 = 1;
# recursion a_{n+1} (n+1)^2 = a_n (n+1/3)(n+2/3)); t_c = Int_0^Del wc/(1 - s) ds (normalised to t_c = Del + ...);
# C1cf = 27 y dt_c/dDel = -wc, so C1i = -1/wc.  Test object: X := rho_c C1i N1B with rho_c = u/(t_c (3 + E u)) (t_c^2 rho_c^2 Delta = const),
# i.e. the analogue of the conifold rho r_{1,1} with the flat derivative C1i y d/dy = d/dt_c.  X = W(u)/t_c, W := -(u N1B)/((3+Eu) wc).
# Conifold analogue (Prop. 6.1): [t^k] sigma_1 in Lam_{k+1}, i.e. NO constant term in z for k >= 0.  Here we print the z = oo constants.
{
my $M = $MU + 2;
my @wcD = (QQ(1)); for my $n (0..$M) { $wcD[$n+1] = $wcD[$n]*(QQ($n)+QQ('1/3'))*(QQ($n)+QQ('2/3'))/(QQ($n+1)**2); }
my @g = map { my $n=$_; my $v=QQ(0); $v += $wcD[$_] for 0..$n; $v } 0..$M;   # wc/(1-Del)
my @tcD = (QQ(0), map { $g[$_-1]/$_ } 1..$M+1);                                 # t_c(Del)
my @Delu = (QQ(0), -$alpha*3, -3*$alpha**2, -$alpha**3, (QQ(0)) x $M);         # Del(u)
@Delu = @Delu[0..$M+1];
sub compose { my ($f,$g,$K)=@_; my @res=(QQ(0))x($K+1); my @p=(QQ(1),(QQ(0))x$K);  # f(g(u)), g(0) = 0
  for my $n (0..$#$f) { last if $n > $K; for my $i (0..$K) { $res[$i] += $f->[$n]*$p[$i]; } my $np = sermul(\@p,$g,$K); @p = @$np; } return \@res; }
my $tcu = compose(\@tcD, \@Delu, $M+1);     # t_c(u) = -E u + ...
my $wcu = compose(\@wcD, \@Delu, $M+1);     # wc(Del(u))
print "\nt_c(u) = ", join(' + ', map { "($tcu->[$_]) u^$_" } 1..3), " + ...;  wc(u) = ", join(' + ', map { "($wcu->[$_]) u^$_" } 0..2), " + ...\n";
# inverse series u(t_c): u = Sum b_n t^n, solve tcu(u(t)) = t
my @ui = (QQ(0), QQ(1)/$tcu->[1], (QQ(0)) x $M);
for my $n (2..$M+1) { my $c = compose($tcu, \@ui, $n); $ui[$n] = -$c->[$n]/$tcu->[1]; }
my $chk = compose($tcu, \@ui, $M+1); my $okinv = ($chk->[1] == 1); $okinv &&= ($chk->[$_] == 0) for 2..$M+1;
print "u(t_c) inverse series check: ", ($okinv?"ok":"FAILED"), "\n";
# 1/((3+Eu) wc(u)) as u-series
my @inv_wc = (QQ(1)/$wcu->[0]); for my $n (1..$M+1) { my $v=QQ(0); $v += $wcu->[$_]*$inv_wc[$n-$_] for 1..$n; $inv_wc[$n] = -$v/$wcu->[0]; }
my $pre = sermul(\@inv_wc, \@inv3Eu, $M+1);     # 1/((3+Eu) wc)
# W(u) and Wr(u) (identity row: rho_c r = u r/(t_c (3+Eu)), C1i absent) as arrays [u^n] -> z-vector
my (@W, @Wr);
for my $n (0..$MU+1) { my @zv = map { my $k=$_; my $v=QQ(0); for my $i (0..$n) { $v += $pre->[$i]*N1f($k,$n-1-$i); } -$v } 0..$N; $W[$n] = zmul(\@zv,\@hinv);
                       my @zr = map { my $k=$_; my $v=QQ(0); for my $i (1..$n) { $v += $inv3Eu[$i-1]*$r[$k][$n-$i]; } $v } 0..$N; $Wr[$n] = zmul(\@zr,\@hinv); }
# compose with u(t_c): [t_c^q] W(u(t_c)) = Sum_n W_n [t^q] u(t)^n
my @pw = ([QQ(1),(QQ(0))x($MU+1)]); for my $n (1..$MU+1) { $pw[$n] = sermul($pw[$n-1], \@ui, $MU+1); }
print "\nH-row in t_c with rho_c and C1i = -1/wc (conifold polarisation), X = rho_c C1i N1B:\n";
for my $m (-1..$MU) { my $q = $m+1; my @a = (QQ(0))x($N+1); for my $n (0..$q) { my $c = $pw[$n][$q]; next if $c == 0; $a[$_] += $c*$W[$n][$_] for 0..$N; }
  if ($m == -1) { print "  [t_c^-1] X: ", join(', ', @a), "\n"; next; } fmt("X (t_c)", \@a, $m, "t_c");
  if ($m == 0) { my @nz = map { -3*$E*$_ } @a; my @p = ($nz[0]*$E, map { $nz[$_]*$E + $nz[$_-1] } 1..$N);
    print "    normalised rho = -3E rho_c (rho(0) = 1): [t_c^-1] rho C1i N1B = E (cf. conifold c/t); [t_c^0] rho C1i N1B * (E+z), z-coeffs: ", join(', ', @p), "\n"; } }
print "\nIdentity row in t_c, rho_c rB:\n";
for my $m (0..$MU) { my $q = $m+1; my @a = (QQ(0))x($N+1); for my $n (0..$q) { my $c = $pw[$n][$q]; next if $c == 0; $a[$_] += $c*$Wr[$n][$_] for 0..$N; }
  fmt("rho_c rB (t_c)", \@a, $m, "t_c"); }
}
# ---- optional dump (env DUMP=1): numerators of [u^m] rB (frame B, point normalisation r = P0/cPP_0), m = 0..MU, for frob_vs_row.pl ----
if ($ENV{DUMP}) { open my $oh, '>', "rB_num_${E1}_${E2}.txt" or die;
  for my $m (0..$MU) { my @c = @{ zmul(zser(sub{ $r[$_[0]][$_[1]] },$m), \@hinv) }; my $J = $m+1;
    for my $j (1..$J) { my @nc = (QQ(0))x($N+1); for my $i (0..$N) { $nc[$i] += $E*$c[$i]; $nc[$i+1] += $j*$c[$i] if $i < $N; } @c = @nc; }
    for my $i ($J+1..$N) { die "dump: not a polynomial at m=$m" if $c[$i] != 0; }
    print $oh "$m ", join(' ', @c[0..$J]), "\n"; }
  close $oh; print "dumped rB numerators to rB_num_${E1}_${E2}.txt (line: m c_0..c_{m+1}; [u^m] rB = Sum c_i z^i / prod_{j<=m+1}(E+jz))\n"; }
