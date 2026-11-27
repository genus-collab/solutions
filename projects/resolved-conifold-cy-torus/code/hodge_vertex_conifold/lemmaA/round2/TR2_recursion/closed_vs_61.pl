# closed_vs_61.pl: compare the closed one-point formula (from ABDKS Prop. 3.1, n = 1)
#   [x^j] F_{g,1} = (1/j^2) [hbar^(2g)] [w^(j+1)] (1+w)^(j(tau+1)) exp( Sum_k hbar^(2k) w^(2k) Q_k(w) )
# with the conjectural formula (6.1) for g = 2..G, j = 1..JMAX, at a given tau (exact BigRat).  usage: perl closed_vs_61.pl tau G JMAX
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $tau = R(shift // 2); my $G = shift // 8; my $JM = shift // 0; my $sig = $tau*($tau+1);
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
sub binomr { my ($a,$i)=@_; return R(0) if $i<0; my $r=R(1); for my $k (0..$i-1){ $r *= ($a-$k); } $r/fact($i) }
my $K=$G+1; my @s = map { R(1)/(R(4)**$_ * fact(2*$_+1)) } 0..$K;
my @sp=(R(1)); for my $n (1..$K){ my $acc=R(0); $acc += $s[$_]*$sp[$n-$_] for 1..$n; push @sp, -$acc; }
my @aa=(R(0),R(-1)); my @al=(R(1)/$tau, R(-1)/($tau+1));
sub ck { my ($u,$a,$k)=@_; my $c=R(0); for my $p (0..$k){ $c += $s[$p]*$sp[$k-$p]*($u**(2*$p))*($a**(2*($k-$p))); } $c }
sub xi_closed { my ($g,$j)=@_; my $W = $j+2;
  my @Q; for my $k (1..$g){ my @q = map { R(0) } 0..$W; for my $i (0..1){ my $cf = fact(2*$k-1)*($j/$al[$i])*ck(R($j),$al[$i],$k);
      for my $l (0..$W){ $q[$l] += $cf*binomr(R(-2*$k),$l)*((-$aa[$i])**$l); } }
    my @qs = map { R(0) } 0..$W; for my $l (0..$W-2*$k){ $qs[$l+2*$k] = $q[$l]; } $Q[$k]=\@qs; }
  my $mul = sub { my ($A,$B)=@_; my @c=map { R(0) } 0..$W; for my $i (0..$W){ next if $A->[$i]->is_zero; for my $l (0..$W-$i){ next if $B->[$l]->is_zero; $c[$i+$l] += $A->[$i]*$B->[$l]; } } \@c };
  my @E=([R(1), map { R(0) } 1..$W]);
  for my $m (1..$g){ my @acc=map { R(0) } 0..$W; for my $k (1..$m){ my $t=$mul->($Q[$k],$E[$m-$k]); $acc[$_] += (R($k)/$m)*$t->[$_] for 0..$W; } $E[$m]=\@acc; }
  my @onew = map { binomr($j*($tau+1),$_) } 0..$W; my $P = $mul->(\@onew,$E[$g]);
  ($P->[$j+1]/($j*$j))*($sig**($g-1))/binomr($j*($tau+1)-1,$j-1) }
my @B=(R(1)); for my $n (1..2*$G+2){ my $acc=R(0); for my $k (0..$n-1){ $acc += binomr(R($n+1),$k)*$B[$k]; } push @B, -$acc/($n+1); }
my @b=(R(0)); for my $k (1..$G+1){ push @b, -$B[2*$k]/(2*$k*fact(2*$k)); }
sub xi61 { my ($h,$j)=@_; my $e1=-$tau; my $e2=$tau+1; my @lg = map { R(0) } 0..$h;
  for my $k (1..$h){ $lg[$k] = $b[$k]*((1+$j*$e2)*($e1**(2*$k)) + (1+$j*$e1)*($e2**(2*$k))); }
  my @ex=(R(1)); for my $m (1..$h){ my $acc=R(0); $acc += $_*$lg[$_]*$ex[$m-$_] for 1..$m; push @ex, $acc/$m; }
  my $pr=R(1); $pr *= ($j-$_) for 1..2*$h-2; $pr*$ex[$h]/($sig**$h) }
my ($n,$bad,$deg)=(0,0,0);
for my $g (2..$G){ my $jm = $JM || 3*$g+1; my @nz; for my $j (1..$jm){
    # skip the special framings where the P2 normalisation binom(j(tau+1)-1, j-1) vanishes (Xi is then not [x^j]F / binom)
    if (binomr($j*($tau+1)-1,$j-1)->is_zero) { $deg++; next; }
    my $a=xi_closed($g,$j); my $c=xi61($g,$j); $n++; if ($a != $c){ $bad++; print "MISMATCH g=$g j=$j closed=$a (6.1)=$c\n"; } push @nz, $j unless $a->is_zero; }
  print "tau=$tau g=$g: j=1..$jm compared; nonzero for j in {".join(",",@nz)."}\n"; }
print "tau=$tau: $n comparisons, $bad mismatches ($deg points with vanishing P2 normalisation skipped)\n";
