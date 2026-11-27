use strict; use warnings; use Math::BigRat;
# frob_vs_row.pl -- exact identification of the frame-B identity row of K_P2 x C^2 at the singular idempotent with the vanishing-cycle
# Frobenius solution of the order-3 operator of frob.pl, for lam = +E and lam = -E (fixes the sign of the leaf scale).
# Conifold model (leaf_structure_proof.md, KEY FORMULA): rho r_{0,1} = (c/(c+z)) (tau/t) 2F1(..; c/z + 2; tau) Exp[-dphi(tau)/z], i.e.
#   log( r (c+z)/c ) - log(hypergeometric series) = alpha(t) + beta(t)/z  (z-dependence only through 1 and 1/z at every t-order).
# Test here: G(u, z) := rB(u, z) (E+z)/E / F(w(u), z),  F = Sum_n b_n(z) w^n the Frobenius series of exponent rho0 (w^rho0 stripped),
#   w = Y - 1 = -(1 + 27 y) = alpha u (3 + 3 alpha u + alpha^2 u^2), alpha = E/3, u = 1/L (as in hrow.pl).
#   For every u-order m = 1..MU check that z [u^m] log G(u, z) is affine in z (alpha_m z + beta_m), exactly, at >= 5 rational z values.
# Input: rB_num_<e1>_<e2>.txt from  DUMP=1 perl hrow.pl cPP_z9.m e1 e2 6  (exact numerators of [u^m] rB over prod_{j<=m+1}(E+jz)).
# Expected: PASS only for lam = +E with the exponent-s solution (s = 1 + E/z), i.e. the sign of frob.pl as written.
# Usage: perl frob_vs_row.pl e1 e2   (output: frob_vs_row_out_<e1>_<e2>.txt by redirection)
my ($E1,$E2) = (shift // 7, shift // 3);
sub QQ { Math::BigRat->new(shift) }
my $E = QQ($E1)+QQ($E2); my $alpha = $E/3;
open my $fh, '<', "rB_num_${E1}_${E2}.txt" or die "run DUMP=1 perl hrow.pl first"; my @num; my $MU = -1;
while (<$fh>) { my @f = split; my $m = shift @f; $num[$m] = [ map { QQ($_) } @f ]; $MU = $m if $m > $MU; } close $fh;
print "eps = ($E1,$E2), E = $E, u-orders 1..$MU\n";
sub padd { my ($a,$b)=@_; my $n = $#$a > $#$b ? $#$a : $#$b; [ map { ($a->[$_]//QQ(0)) + ($b->[$_]//QQ(0)) } 0..$n ] }
sub psc  { my ($a,$c)=@_; [ map { $_*$c } @$a ] }
sub pmul { my ($a,$b)=@_; my @c = (QQ(0)) x (@$a+@$b-1); for my $i (0..$#$a) { next if $a->[$i]==0; for my $j (0..$#$b) { $c[$i+$j] += $a->[$i]*$b->[$j]; } } \@c }
sub peval { my ($a,$x)=@_; my $v=QQ(0); for my $i (reverse 0..$#$a) { $v = $v*$x + $a->[$i]; } $v }
sub bseries { my ($lam,$which,$NMAX)=@_;   # returns ([N_n], [D_n]) polynomials in zeta = 1/z, b_n = N_n/D_n
  my $lamz = [QQ(0), QQ($lam)]; my @a = map { padd(psc($lamz, QQ(-1)/3), [QQ($_)/3]) } (0,1,2);
  my $rho0 = $which eq 's' ? padd([QQ(1)], $lamz) : [QQ(1)];
  my $act = sub { my ($rho)=@_;
    my $th = sub { my ($f)=@_; my %g; for my $k (keys %$f) { my $e = padd($rho,[QQ($k)]); my $t = pmul($e,$f->{$k});
        $g{$k-1} = padd($g{$k-1}//[QQ(0)], $t); $g{$k} = padd($g{$k}//[QQ(0)], $t); } \%g };
    my $t3 = $th->($th->($th->({0=>[QQ(1)]})));
    my $p = {0=>[QQ(1)]}; for my $ai (@a) { my $q = $th->($p); for my $k (keys %$p) { $q->{$k} = padd($q->{$k}//[QQ(0)], pmul($ai,$p->{$k})); } $p = $q; }
    my %r = %$t3; for my $k (keys %$p) { $r{$k} = padd($r{$k}//[QQ(0)], psc($p->{$k},QQ(-1))); $r{$k+1} = padd($r{$k+1}//[QQ(0)], psc($p->{$k},QQ(-1))); }
    \%r };
  my %Cc; my $C = sub { my ($n,$j)=@_; $Cc{"$n,$j"} //= do { my $r = $act->(padd($rho0,[QQ($n)])); $r->{$j-2} // [QQ(0)] }; };
  my @Nn = ([QQ(1)]); my @Dn = ([QQ(1)]);
  for my $n (1..$NMAX) { my $acc = [QQ(0)];
    for my $j (1..3) { next if $n-$j < 0; my $t = pmul($C->($n-$j,$j), $Nn[$n-$j]); for my $m ($n-$j+1..$n-1) { $t = pmul($t, $C->($m,0)); } $acc = padd($acc,$t); }
    $Nn[$n] = psc($acc, QQ(-1)); $Dn[$n] = pmul($Dn[$n-1], $C->($n,0)); }
  return (\@Nn, \@Dn); }
sub sermul { my ($a,$b,$M)=@_; my @c=(QQ(0))x($M+1); for my $i(0..$M){for my $j(0..$M-$i){$c[$i+$j]+=($a->[$i]//QQ(0))*($b->[$j]//QQ(0))}} return \@c; }
sub serlog { my ($g,$M)=@_; # g_0 = 1; log g via (log g)' = g'/g
  my @inv = (QQ(1)); for my $n (1..$M) { my $v=QQ(0); $v += $g->[$_]*$inv[$n-$_] for 1..$n; $inv[$n] = -$v; }
  my @dg = map { ($_+1)*$g->[$_+1] } 0..$M-1; my $q = sermul(\@dg, \@inv, $M-1); return [QQ(0), map { $q->[$_-1]/$_ } 1..$M]; }
my @w = (QQ(0), $alpha*3, 3*$alpha**2, $alpha**3); push @w, (QQ(0)) x $MU; @w = @w[0..$MU];
my @wp = ([QQ(1),(QQ(0))x$MU]); for my $n (1..$MU) { $wp[$n] = sermul($wp[$n-1], \@w, $MU); }
my @zc = grep { my $z = $_; my $bad = 0; for my $j (1..12) { $bad = 1 if $z*$j == $E || $z*$j == -$E; } !$bad }
         map { QQ($_) } ('1/3','2/7','7/3','13/4','11','17/5','19/6','23/9','29/11','31/7','37/8','41/13');
print "z sample points: @zc\n";
for my $sg (1,-1) { for my $which ('s','1') {
  my $lam = $sg*$E; my ($Nn,$Dn) = bseries($lam,$which,$MU);
  my %ell;   # ell{m} = list of z*[u^m] log G at the sample points
  for my $z0 (@zc) { my $zeta = QQ(1)/$z0;
    my @b = map { peval($Nn->[$_],$zeta)/peval($Dn->[$_],$zeta) } 0..$MU;
    my @F = (QQ(0))x($MU+1); for my $n (0..$MU) { $F[$_] += $b[$n]*$wp[$n][$_] for 0..$MU; }
    my @R = map { my $m=$_; my $d=QQ(1); $d *= ($E + $_*$z0) for 1..$m+1; peval($num[$m],$z0)/$d*($E+$z0)/$E } 0..$MU;
    die "R_0 != 1" unless $R[0] == 1 && $F[0] == 1;
    my @Finv = (QQ(1)); for my $n (1..$MU) { my $v=QQ(0); $v += $F[$_]*$Finv[$n-$_] for 1..$n; $Finv[$n] = -$v; }
    my $G = sermul(\@R, \@Finv, $MU); my $lg = serlog($G, $MU);
    push @{$ell{$_}}, $z0*$lg->[$_] for 1..$MU; }
  my @verdict;
  for my $m (1..$MU) { my @v = @{$ell{$m}}; my $sl = ($v[1]-$v[0])/($zc[1]-$zc[0]); my $ok = 1;
    for my $i (2..$#zc) { $ok = 0 unless $v[$i] == $v[0] + $sl*($zc[$i]-$zc[0]); }
    push @verdict, "u^$m:".($ok ? "affine (alpha=".$sl.", beta=".($v[0]-$sl*$zc[0]).")" : "NOT affine"); }
  print "lam = $lam, solution w^", ($which eq 's' ? 's' : '1'), ":\n  ", join("\n  ", @verdict), "\n";
} }
