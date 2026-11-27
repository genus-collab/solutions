# onepoint.pl: the one-point function of C_tau in closed form (consequence of ABDKS Prop. 3.1, n = 1, after integrating by parts):
#   [x^j] F_{g,1} = -(1/j^2) [hbar^(2g)] Res_{z=oo} exp(-j S(j hbar d_z) xhat(z)) dz,
#   -j S(j hbar d) xhat = - Sum_a (j/alpha_a) (S(j hbar d)/S(alpha_a hbar d)) log(z - a),  (a, 1/alpha_a) = (0, tau), (-1, -(tau+1)).
# Expanded at z = oo (w = 1/z):  exp(...) = w^(-j) (1+w)^(j(tau+1)) exp( Sum_{k>=1} hbar^(2k) w^(2k) Q_k(w) ),
#   Q_k(w) = (2k-1)! Sum_a (j/alpha_a) c_k(j,alpha_a) (1 - a w)^(-2k),  c_k(u,al) = [t^(2k)] S(u t)/S(al t),
# and  [x^j] F_{g,1} = (1/j^2) [hbar^(2g)] [w^(j+1)] (1+w)^(j(tau+1)) exp(...)   (only w^(>=2g) at hbar^(2g): zero for j <= 2g-2).
# Also evaluates (6.1): sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j-l) [lam^(2h)] S1(e1 lam)^(1+j e2) S1(e2 lam)^(1+j e1), S1 = (y/2)/sinh(y/2).
# Compares with the full formula (XY.pm, which agrees with exact TR).  usage: perl onepoint.pl tau G JX
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use XY qw(xi_formula_n1);
sub R { Math::BigRat->new(@_) }
my $tau = R(shift // 2); my $G = shift // 5; my $JX = shift // 3; my $sig = $tau*($tau+1);
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
sub binomr { my ($a,$i)=@_; return R(0) if $i<0; my $r=R(1); for my $k (0..$i-1){ $r *= ($a-$k); } $r/fact($i) }
my $K=$G+1; my @s = map { R(1)/(R(4)**$_ * fact(2*$_+1)) } 0..$K;            # S(t) = 2 sinh(t/2)/t = Sum s_k t^2k
my @sp=(R(1)); for my $n (1..$K){ my $acc=R(0); $acc += $s[$_]*$sp[$n-$_] for 1..$n; push @sp, -$acc; }   # 1/S
my @aa=(R(0),R(-1)); my @al=(R(1)/$tau, R(-1)/($tau+1));
sub ck { my ($u,$a,$k)=@_; my $c=R(0); for my $p (0..$k){ $c += $s[$p]*$sp[$k-$p]*($u**(2*$p))*($a**(2*($k-$p))); } $c }
sub xi_closed { my ($g,$j)=@_; my $W = $j+2; # w-series length
  # series in (h2 = hbar^2 power m, w power): exp( Sum_k h2^k w^(2k) Q_k(w) )
  my @Q; for my $k (1..$g){ my @q = map { R(0) } 0..$W; for my $i (0..1){ my $cf = fact(2*$k-1)*($j/$al[$i])*ck(R($j),$al[$i],$k);
      for my $l (0..$W){ $q[$l] += $cf*binomr(R(-2*$k),$l)*((-$aa[$i])**$l); } }   # (1 - a w)^(-2k)
    my @qs = map { R(0) } 0..$W; for my $l (0..$W-2*$k){ $qs[$l+2*$k] = $q[$l]; } $Q[$k]=\@qs; }
  my $mul = sub { my ($A,$B)=@_; my @c=map { R(0) } 0..$W; for my $i (0..$W){ next if $A->[$i]->is_zero; for my $l (0..$W-$i){ next if $B->[$l]->is_zero; $c[$i+$l] += $A->[$i]*$B->[$l]; } } \@c };
  my @E=([R(1), map { R(0) } 1..$W]);
  for my $m (1..$g){ my @acc=map { R(0) } 0..$W; for my $k (1..$m){ my $t=$mul->($Q[$k],$E[$m-$k]); $acc[$_] += (R($k)/$m)*$t->[$_] for 0..$W; } $E[$m]=\@acc; }
  my @onew = map { binomr($j*($tau+1),$_) } 0..$W; my $P = $mul->(\@onew,$E[$g]);
  my $c = $P->[$j+1]/($j*$j);       # [x^j] F_{g,1}
  $c*($sig**($g-1))/binomr($j*($tau+1)-1,$j-1) }
# (6.1)
my @b; { my @B=(R(1)); for my $n (1..2*$G+2){ my $acc=R(0); for my $k (0..$n-1){ $acc += binomr(R($n+1),$k)*$B[$k]; } push @B, -$acc/($n+1); }
  # log S1(y) = log((y/2)/sinh(y/2)) = -Sum_{k>=1} B_{2k} y^{2k}/(2k (2k)!)
  @b = (R(0)); for my $k (1..$G+1){ push @b, -$B[2*$k]/(2*$k*fact(2*$k)); } }
sub xi61 { my ($h,$j)=@_; my $e1=-$tau; my $e2=$tau+1; my @lg = map { R(0) } 0..$h; # log of S1(e1 l)^(1+j e2) S1(e2 l)^(1+j e1), in lam^2
  for my $k (1..$h){ $lg[$k] = $b[$k]*((1+$j*$e2)*($e1**(2*$k)) + (1+$j*$e1)*($e2**(2*$k))); }
  my @ex=(R(1)); for my $m (1..$h){ my $acc=R(0); $acc += $_*$lg[$_]*$ex[$m-$_] for 1..$m; push @ex, $acc/$m; }
  my $pr=R(1); $pr *= ($j-$_) for 1..2*$h-2; $pr*$ex[$h]/($sig**$h) }
my $XF = xi_formula_n1($tau,$G,2*$G-2+$JX); my ($n,$bad1,$bad2)=(0,0,0);
for my $g (1..$G){ my @row; for my $j (1..2*$g-2+$JX){ my $a=xi_closed($g,$j); my $f=$XF->{"$g,$j"}; my $c = $g>=2 ? xi61($g,$j) : undef; $n++;
    $bad1++ unless $a==$f; $bad2++ if defined $c && $c != $a;
    push @row, "$j:".($a->is_zero ? "0" : $a).($a==$f ? "" : "[!=full $f]").(defined $c && $c!=$a ? "[!=(6.1) $c]" : ""); }
  print "tau=$tau g=$g: ", join("  ",@row), "\n"; }
print "tau=$tau: $n values; closed vs full formula: $bad1 mismatches; closed vs (6.1) (g>=2): $bad2 mismatches\n";
