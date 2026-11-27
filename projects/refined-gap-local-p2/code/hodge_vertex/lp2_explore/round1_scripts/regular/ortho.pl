use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RegLib; use Math::BigRat;
# Test 7, part C: rows at all three colours (incl. the two regular ones at the conifold), exact in Q(omega), and the
# orthogonality relations Pi_ij(z) := Sum_ab P_a^i(z) eta^ab P_b^j(-z) = delta_ij  mod z^(K+1).
# Rows exactly as lp2_lib.m buildRmatrix: P0 = Sum cPP_k z^k (L -> L(omega^i t)), P1 = (L P0 + z th P0) C1i,
# P2 = (L P1 + z t^3 Dy P1 - d2 P1)/c2 with c2 = C1i^2/(1+27y), d2 = E (C1i^2 - 1)/(2(1+27y)); t^3 Dy P1 = C1i th Q1 - XX C1i Q1.
# C1i and XX are set to arbitrary rationals (two choices) to show they drop out of Pi.
# Usage: perl ortho.pl e1 e2 t0 [C1i XX]  (t0 rational; t0 = -1/3 is the conifold on the real branch: colour 0 is skipped there)
my ($E1,$E2,$T0,$C1I,$XXV) = @ARGV; $E1=QQ($E1); $E2=QQ($E2); my $E=$E1+$E2; my $al=$E/3; my $t0 = QQ($T0);
$C1I = QQ($C1I // '3/7'); $XXV = QQ($XXV // '-2/9');
my @S; open my $in, '<', "$FindBin::Bin/Sk_${E1}_${E2}.txt" or die; while (<$in>) { chomp; my ($k,@t) = split ' '; my %p; for (@t) { my ($e,$c) = split /:/; $p{$e} = QQ($c) } $S[$k] = \%p } close $in;
my $K = $#S;
if ($ENV{PERT}) { my ($c,$j,$d) = split /,/, $ENV{PERT}; $S[$c] = lp_add($S[$c], {$j => QQ($d)}); print "CONTROL: S_$c += $d s^$j
"; }
# eta, its inverse (exact 3x3)
my $e3 = $al**3; my $kk = 1/(3*$E1*$E2);
my @eta = ([ 1/(3*$e3*$E1*$E2), 0, 0 ], [ 0, 0, -$kk ], [ 0, -$kk, $E*$kk ]);   # etaF = -{{-1/(3 e3 e1 e2),0,0},{0,0,k},{0,k,-E k}}
$_ = [ map { QQ($_) } @$_ ] for @eta;
sub inv3 { my $m = shift; my @a = map { [@$_] } @$m;
  my $det = $a[0][0]*($a[1][1]*$a[2][2]-$a[1][2]*$a[2][1]) - $a[0][1]*($a[1][0]*$a[2][2]-$a[1][2]*$a[2][0]) + $a[0][2]*($a[1][0]*$a[2][1]-$a[1][1]*$a[2][0]);
  my @c; for my $i (0..2) { for my $j (0..2) { my @r = grep { $_ != $j } 0..2; my @cc = grep { $_ != $i } 0..2;
     my $minor = $a[$r[0]][$cc[0]]*$a[$r[1]][$cc[1]] - $a[$r[0]][$cc[1]]*$a[$r[1]][$cc[0]];
     $c[$i][$j] = (($i+$j)%2 ? -1 : 1) * $minor / $det; } } return \@c; }
my $ei = inv3(\@eta);
print "eta^-1 = [", join('; ', map { join(', ', @$_) } @$ei), "]   (expect eta^00 = 3 e1 e2 alpha^3 = ", 3*$E1*$E2*$e3, ", eta^11 = -3 e1 e2 E, eta^12 = -3 e1 e2, eta^22 = 0)\n";
# z-series of Laurent polys -> Q(omega) at s
my $Y = {-1 => QQ(1), 0 => -$al};
sub thS { my $F = shift; [ map { lp_theta($_,$al) } @$F ] }
my $th0 = thS(\@S);
my @Q1 = map { my $m=$_; my $r = lp_mul($Y,$S[$m]); $r = lp_add($r, $th0->[$m-1]) if $m>=1; $r } 0..$K;    # Q1 = N P0
my $thQ = thS(\@Q1);
my @NQ = map { my $m=$_; my $r = lp_mul($Y,$Q1[$m]); $r = lp_add($r, $thQ->[$m-1]) if $m>=1; $r } 0..$K;   # N Q1
my $w = [QQ(0),QQ(1)]; my $one = [QQ(1),QQ(0)];
my $pref = 3*$E**3/($E1*$E2);   # Sqrt3^2 sq^-2
my $oneplus27y = 1 + 27*$t0**3;
my %rows; my @cols = ($t0 == QQ(-1)/3) ? (1,2) : (0,1,2);
for my $i (@cols) {
  my $tau = w_scale(w_pow($w,$i), $t0);                       # omega^i t0
  my $num = w_add($one, w_scale($tau,QQ(3)));                  # 1 + 3 tau
  my $si = w_mul($num, w_inv(w_scale($tau, 3*$al)));          # s = (1+3tau)/(3 alpha tau)
  my $qi = w_sub($one, w_scale($si,$al)); my $chk = w_sub($one, w_inv(w_pow($qi,3)));   # 1+27y = 1 - q^-3
  die "1+27y mismatch" unless w_eq($chk, w_fromq($oneplus27y));
  my $Li = w_sub(w_inv($si), w_fromq($al));
  my @p0 = map { lp_eval($S[$_],$si) } 0..$K; my @q1 = map { lp_eval($Q1[$_],$si) } 0..$K; my @nq = map { lp_eval($NQ[$_],$si) } 0..$K;
  my @p1 = map { w_scale($_, $C1I) } @q1;
  my @p2; for my $m (0..$K) { my $v = w_scale($nq[$m], $oneplus27y); $v = w_add($v, w_scale($q1[$m], $E/2));
     $v = w_sub($v, w_scale($q1[$m-1], $oneplus27y*$XXV)) if $m >= 1;
     $p2[$m] = w_sub(w_scale($v, 1/$C1I), w_scale($q1[$m], $C1I*$E/2)); }
  $rows{$i} = [\@p0,\@p1,\@p2];
  printf "colour %d: s = %s, L = %s, L+alpha = %s; |P_a[m]| finite for all m <= %d (exact Q(omega) values)\n", $i, w_str($si), w_str($Li), w_str(w_add($Li,w_fromq($al))), $K;
  print "   P0[0..2] = ", join(' | ', map { w_str(w_scale($p0[$_],1)) } 0..2), "  (times Sqrt3/sq)\n";
}
my $allok = 1;
for my $i (@cols) { for my $j (@cols) {
  my @Pi;
  for my $m (0..$K) { my $acc = [QQ(0),QQ(0)];
    for my $n (0..$m) { my $sg = ($m-$n) % 2 ? -1 : 1;   # P_b^j(-z): z^(m-n) coefficient times (-1)^(m-n)
      for my $a (0..2) { for my $b (0..2) { next if $ei->[$a][$b] == 0;
        $acc = w_add($acc, w_scale(w_mul($rows{$i}[$a][$n], $rows{$j}[$b][$m-$n]), $ei->[$a][$b]*$sg*$pref)); } } }
    push @Pi, $acc; }
  my $ok = 1; for my $m (0..$K) { my $want = ($i==$j && $m==0) ? $one : [QQ(0),QQ(0)]; $ok = 0 unless w_eq($Pi[$m],$want) }
  $allok &&= $ok;
  print "Pi_{$i$j}: ", ($ok ? "= delta_ij exactly through z^$K" : "DEVIATES: ".join(' ; ', map { "z^$_: ".w_str($Pi[$_]) } 0..$K)), "\n";
} }
print "ALL ORTHOGONALITY RELATIONS: ", ($allok ? "EXACT" : "FAIL"), "  (t0 = $t0, C1i = $C1I, XX = $XXV)\n";
