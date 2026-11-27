use strict; use warnings; use Math::BigRat;
# Leaf-support test for the frame-B identity row of K_P2 x C^2 at the singular idempotent near the conifold.
# cPP_k (file cPP_z9.m) = Sqrt[3]^(+-1) sq^-1 N_k(L)/(E+3L)^(k+1) with deg N_k <= k (checked below).
# Normalised row r_k(L) = cPP_k/cPP_0; frame B: rowB(z) = Sum_k r_k z^k * h_w(z)^(-1), w = (E,-e1,-e2) (float test: = E/(E+z) at L = oo).
# Expand in u = 1/L (u = -(1+3t)/alpha, a unit times t_c), and test whether [u^m] rowB(z) * Prod_{j=1}^{m+1} (E + j z) is a
# polynomial in z of degree <= m+1 (printing the coefficients z^0..z^9).
my $f = shift; my ($E1,$E2) = (shift // 7, shift // 3); my $MU = shift // 4;
local $/; open my $fh,'<',$f or die; my $s=<$fh>; close $fh;
$s =~ s/\s+//g; $s =~ s/^\{//; $s =~ s/\}$//;
my @parts; my $d=0; my $cur='';
for my $ch (split //,$s){ if($ch eq '['||$ch eq '('){$d++} elsif($ch eq ']'||$ch eq ')'){$d--}
  if($ch eq ',' && $d==0){push @parts,$cur;$cur=''} else {$cur.=$ch} } push @parts,$cur;
my $N = $#parts;
my $bs = chr(92);
our ($e1,$e2,$s3,$sqv,$L);
sub QQ { Math::BigRat->new(shift) }
$e1 = QQ($E1); $e2 = QQ($E2); my $E = $e1+$e2;
# exact evaluation with Sqrt[3] -> s3 = 1 and sq -> 1; the true value differs by 3^(p/2) * sq^-1, p = +-1; fixed by the float check
sub ev { my ($x,$Lv,$s3v)=@_; my $t=$x; $t =~ s/\Q$bs\E//g;
  $t =~ s/\Q[Epsilon][1]\E/\$e1/g; $t =~ s/\Q[Epsilon][2]\E/\$e2/g; $t =~ s/\QLLL[0][y]\E/\$L/g;
  $t =~ s/\QSqrt[3]\E/\$s3/g; $t =~ s/\bsq\b/\$sqv/g; $t =~ s/\^/**/g; $t =~ tr/[]/()/;
  # make integer literals BigRat
  $t =~ s/(?<![\$\w.])(\d+)/QQ('$1')/g;
  local $L = $Lv; local $s3 = $s3v; local $sqv = QQ(1);
  my $v = eval $t; die "$@ in $t" if $@; return $v; }
# float version for the sqrt(3) power: compare ratio with s3 = 1 vs s3 = sqrt3 numerically
my @fac;
for my $k (0..$N) {
  my $a = ev($parts[$k], QQ(5), QQ(1)); my $b = ev($parts[$k], QQ(5), QQ(3));   # s3 = 3 gives 3^p
  my $q = $b/$a; # = 3^p
  my $p = ($q == 3) ? 1 : ($q == QQ(1)/3 ? -1 : ($q==1?0:'?')); push @fac, $p; }
print "Sqrt[3] powers: @fac\n";
# exact N_k(L) = cPP_k (E+3L)^(k+1) (with s3 = 1, sq = 1), degree check by interpolation at k+3 points
my @Ncoef;
for my $k (0..$N) {
  my @xs = map { QQ($_+2) } 0..($k+2); my @ys = map { ev($parts[$k], $_, QQ(1)) * ($E + 3*$_)**($k+1) } @xs;
  # Newton divided differences -> monomial coefficients
  my @c = @ys; my $n = scalar(@xs);
  for my $j (1..$n-1) { for (my $i=$n-1; $i>=$j; $i--) { $c[$i] = ($c[$i]-$c[$i-1])/($xs[$i]-$xs[$i-$j]); } }
  my @poly = (QQ(0)) x $n; # convert Newton form
  my @basis = (QQ(1));
  for my $i (0..$n-1) { for my $m (0..$#basis) { $poly[$m] += $c[$i]*$basis[$m]; }
    my @nb = (QQ(0)) x (scalar(@basis)+1); for my $m (0..$#basis) { $nb[$m+1] += $basis[$m]; $nb[$m] -= $xs[$i]*$basis[$m]; } @basis = @nb; }
  my $deg = 0; for my $m (0..$n-1) { $deg = $m if $poly[$m] != 0; }
  print "k=$k deg N_k = $deg\n" if $deg > $k;
  $Ncoef[$k] = [ @poly[0..$k] ];
}
# r_k(u) = 3^((p_k - p_0)/2) [Sum_i n_{k,i} L^i]/(n_0 (E+3L)^k) = 3^.. [Sum_i n_{k,i} u^(k-i)] / (n_0 (3 + E u)^k)
# series in u to order MU
sub sermul { my ($a,$b)=@_; my @c=(QQ(0))x($MU+1); for my $i(0..$MU){for my $j(0..$MU-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}} return \@c; }
my @inv3Eu = map { (QQ(-1)*$E)**$_ / QQ(3)**($_+1) } 0..$MU;   # 1/(3+Eu)
my @rowu; # rowu[k] = arrayref of u-coefficients
for my $k (0..$N) {
  my $pp = $fac[$k] - $fac[0]; die "odd sqrt3 power" if $pp % 2; my $pref = QQ(3)**($pp/2);
  my @num = (QQ(0))x($MU+1); for my $i (0..$k) { my $e = $k-$i; $num[$e] = $Ncoef[$k][$i] if $e <= $MU; }
  my $ser = \@num; for (1..$k) { $ser = sermul($ser, \@inv3Eu); }
  $rowu[$k] = [ map { $_ * $pref / $Ncoef[0][0] } @$ser ];
}
# h_w^(-1): log h = -Sum B2k/(2k(2k-1)) z^(2k-1) S_{2k-1}
my @B = (QQ('1/6'),QQ('-1/30'),QQ('1/42'),QQ('-1/30'),QQ('5/66'));
my @w = ($E, -$e1, -$e2);
my @logh = (QQ(0)) x ($N+1);
for my $k (1..5) { my $m = 2*$k-1; next if $m > $N; my $S = QQ(0); $S += QQ(1)/($_**$m) for @w; $logh[$m] = -$B[$k-1]/(2*$k*(2*$k-1))*$S; }
my @neg = map { -$_ } @logh;
my @hinv = (QQ(1), (QQ(0)) x $N);
for my $n (1..$N) { my $v = QQ(0); for my $k (1..$n) { $v += $k*$neg[$k]*$hinv[$n-$k]; } $hinv[$n] = $v/$n; }
for my $m (0..$MU) {
  my @a = map { $rowu[$_][$m] } 0..$N;            # z-series of [u^m] row
  my @b = (QQ(0)) x ($N+1); for my $i (0..$N) { for my $j (0..$N-$i) { $b[$i+$j] += $a[$i]*$hinv[$j]; } }
  my @c = @b; for my $j (1..$m+1) { my @nc = (QQ(0))x($N+1); for my $i (0..$N) { $nc[$i] += $E*$c[$i]; $nc[$i+1] += $j*$c[$i] if $i < $N; } @c = @nc; }
  my $maxnz = -1; for my $i (0..$N) { $maxnz = $i if $c[$i] != 0; }
  print "u^$m: [u^m]rowB * prod_{j<=m+1}(E+jz): highest nonzero z-power = $maxnz (of $N); low coeffs: ",
     join(', ', map { $c[$_] } 0..$N), "\n";
}
