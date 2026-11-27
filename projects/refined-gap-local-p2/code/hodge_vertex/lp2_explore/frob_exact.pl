use strict; use warnings; use Math::BigRat;
# frob_exact.pl -- exact (Math::BigRat) version of frob.pl, run for lam = +E and lam = -E, to fix the sign of the leaf scale once
# against the cPP denominators (the frame-B rows of cPP_z9.m have [u^m] in span{1/(E + j z) : j <= m+1}, simple poles at z = -E/j,
# every j <= m+1 present: hrow.pl, leafrow.pl).
# Operator (frob.pl): theta^3 - Y prod_{i=1..3}(theta + a_i), a_i = -lam/(3z) + (i-1)/3, theta = Y d/dY, w = Y - 1.
#   Exponents at w = 0: {0, 1, s}, s = 1 + lam/z.  Frobenius solutions w^rho0 Sum_n b_n w^n for rho0 = 1 and rho0 = s
#   (rho0 = 0 is resonant with 1 and is skipped).  Everything is polynomial in zeta = 1/z; b_n = N_n(zeta)/D_n(zeta) with
#   D_n = prod_{m=1}^n C_0(rho0 + m); the script cancels common linear factors exactly and prints the pole set in z (with orders).
# Expected (if the row is the vanishing-cycle solution w^s, as in the conifold KEY FORMULA): for lam = +E the exponent-s solution has
#   only simple poles at z = -E/j, j <= n+1; lam = -E puts them at z = +E/j (wrong sign).
# Usage: perl frob_exact.pl [lam_abs NMAX]   (default 10 6).  Output: frob_exact_out.txt (redirect).
my $LAM = shift // 10; my $NMAX = shift // 6;
sub QQ { Math::BigRat->new(shift) }
# polynomials in zeta: array refs of BigRat, index = degree
sub padd { my ($a,$b)=@_; my $n = $#$a > $#$b ? $#$a : $#$b; [ map { ($a->[$_]//QQ(0)) + ($b->[$_]//QQ(0)) } 0..$n ] }
sub psc  { my ($a,$c)=@_; [ map { $_*$c } @$a ] }
sub pmul { my ($a,$b)=@_; my @c = (QQ(0)) x (@$a+@$b-1); for my $i (0..$#$a) { next if $a->[$i]==0; for my $j (0..$#$b) { $c[$i+$j] += $a->[$i]*$b->[$j]; } } \@c }
sub ptrim { my ($a)=@_; my @c=@$a; pop @c while @c>1 && $c[-1]==0; \@c }
sub peval { my ($a,$x)=@_; my $v=QQ(0); for my $i (reverse 0..$#$a) { $v = $v*$x + $a->[$i]; } $v }
sub pdivlin { my ($a,$r)=@_; # divide by (zeta - r), assumes exact
  my @c = @$a; my $n = $#c; my @q = (QQ(0)) x $n; my $carry = QQ(0);
  for (my $i=$n; $i>=1; $i--) { $carry = $c[$i] + $carry*$r; $q[$i-1] = $carry; } return \@q; }
sub run { my ($lam,$which)=@_;
  my $lamz = [QQ(0), QQ($lam)];                         # lam*zeta
  my @a = map { padd(psc($lamz, QQ(-1)/3), [QQ($_)/3]) } (0,1,2);
  my $s = padd([QQ(1)], $lamz);                           # s = 1 + lam zeta
  my $rho0 = $which eq 's' ? $s : [QQ(1)];
  # act(rho): coefficients (polys) of w^(rho + k), k = -3..1
  my $act = sub { my ($rho)=@_;
    my $th = sub { my ($f)=@_; my %g; for my $k (keys %$f) { my $e = padd($rho,[QQ($k)]); my $t = pmul($e,$f->{$k});
        $g{$k-1} = padd($g{$k-1}//[QQ(0)], $t); $g{$k} = padd($g{$k}//[QQ(0)], $t); } \%g };
    my $f = {0=>[QQ(1)]}; my $t3 = $th->($th->($th->($f)));
    my $p = {0=>[QQ(1)]}; for my $ai (@a) { my $q = $th->($p); for my $k (keys %$p) { $q->{$k} = padd($q->{$k}//[QQ(0)], pmul($ai,$p->{$k})); } $p = $q; }
    my %r = %$t3; for my $k (keys %$p) { $r{$k} = padd($r{$k}//[QQ(0)], psc($p->{$k},QQ(-1))); $r{$k+1} = padd($r{$k+1}//[QQ(0)], psc($p->{$k},QQ(-1))); }
    return \%r; };
  my $C = sub { my ($n,$j)=@_; my $r = $act->(padd($rho0,[QQ($n)])); return ptrim($r->{$j-2} // [QQ(0)]); };
  # sanity: offset -3 vanishes
  { my $r = $act->(padd($rho0,[QQ(2)])); my $z3 = ptrim($r->{-3}//[QQ(0)]); die "offset -3 nonzero" unless @$z3==1 && $z3->[0]==0; }
  my @C0 = map { $C->($_,0) } 0..$NMAX;
  my @Nn = ([QQ(1)]); my @Dfac;   # D_n = prod_{m=1}^n C0[m]
  print "lam = $lam, solution w^", ($which eq 's' ? "s (s = 1 + lam/z)" : "1"), ":\n";
  for my $n (1..$NMAX) {
    my $acc = [QQ(0)];
    for my $j (1..3) { next if $n-$j < 0; my $t = pmul($C->($n-$j,$j), $Nn[$n-$j]);
      for my $m ($n-$j+1..$n-1) { $t = pmul($t, $C0[$m]); } $acc = padd($acc,$t); }
    $Nn[$n] = ptrim(psc($acc, QQ(-1)));
  }
  # pole analysis: candidate roots of all C0[m]: linear factors in zeta; C0 is cubic in rho, i.e. a polynomial in zeta
  for my $n (1..$NMAX) {
    my $Nr = $Nn[$n]; my %mult;   # root -> multiplicity in D_n
    my @roots;
    for my $m (1..$n) { my $p = $C0[$m];
      # roots of p: candidates zeta = q/lam for integer q in [-2n-4, 2n+4] (all roots are of this form)
      my $pp = $p; for my $q (-3*$NMAX-6 .. 3*$NMAX+6) { my $r = QQ($q)/QQ($lam); while (@$pp > 1 && peval($pp,$r) == 0) { $mult{"$r"}++; $pp = pdivlin($pp,$r); } }
      die "C0[$m] not split over candidates: deg ".$#$pp if $#$pp > 0; }
    my @out;
    for my $rk (sort { Math::BigRat->new($a) <=> Math::BigRat->new($b) } keys %mult) { my $r = QQ($rk); my $k = $mult{$rk}; my $num = $Nr;
      my $cancel = 0; while ($cancel < $k && @$num > 1 && peval($num,$r) == 0) { $num = pdivlin($num,$r); $cancel++; }
      my $ord = $k - $cancel; next if $ord <= 0;
      my $zp = ($r == 0) ? 'oo' : (QQ(1)/$r); push @out, "z=$zp".($ord>1?"(order $ord)":""); }
    # behaviour at zeta = oo (z = 0): deg N - deg D
    my $degD = 0; $degD += $#{ptrim($C0[$_])} for 1..$n;
    printf "  b_%d: poles %s;  deg_zeta N - deg_zeta D = %d\n", $n, (@out ? join(' ', @out) : 'none'), $#$Nr - $degD;
  }
}
for my $sg (1,-1) { for my $which ('s','1') { run($sg*$LAM, $which); } }
