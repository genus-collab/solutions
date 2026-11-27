# OS hook formula for psi = Y^{-1}, y(z) = Y(z) = (e^z-1-z)/z; compare with Xi_g(mu) mu^mu/mu!
use strict; use warnings; use lib '.'; use ser; use Math::BigRat; use Math::BigInt;
sub R { Math::BigRat->new(@_) }
sub fact { my $n=shift; my $f=Math::BigInt->new(1); $f*=$_ for 2..$n; $f }
my $GMAX = 3; my $MUMAX = 8;
my $NS = 2*$GMAX-1+$MUMAX+2;
$ser::N = $NS;
my @s = map { R(1, fact($_+1)) } 0..$NS;  $s[0]=R(0);   # Y = sum_{k>=1} z^k/(k+1)!
# psi = inverse of Y: Yhat = 2Y = z + ..., revert, psi(y) = Rv(2y)
my $Yh = ser::zero(); for my $k (1..$NS){ $Yh->[$k] = 2*$s[$k]; }
my $Rv = ser::revert($Yh);
my @psi = map { $Rv->[$_] * R(2)**$_ } 0..$NS;
print "psi = ", join(", ", map { $psi[$_]->bstr } 0..5), "\n";
sub Xi { my ($g,$m)=@_; my $j=R($m);
  return ($j-2)/24 if $g==1;
  return ($j-1)*($j-2)*($j-4)*(5*$j-6)/5760 if $g==2;
  return ($j-6)*($j-4)*($j-3)*($j-2)*($j-1)*(35*$j*$j-126*$j+80)/2903040 if $g==3; }
for my $mu (1..$MUMAX) {
  my $D = $mu-1;
  # Q[i][j] = s_{i+j+1}
  my @Q; for my $i (0..$D){ for my $j (0..$D-$i){ $Q[$i][$j] = $s[$i+$j+1]; } }
  my @A; my @P; # P = Q^m (u-u')^(m-1)
  my @cur; for my $i (0..$D){ for my $j (0..$D-$i){ $cur[$i][$j] = ($i==0&&$j==0)?R(1):R(0); } }
  my @Am;
  for my $m (1..$mu) {
    # cur <- cur * Q * (m>1 ? (u-u') : 1)
    my @nw; for my $i (0..$D){ for my $j (0..$D-$i){ $nw[$i][$j]=R(0); } }
    for my $i (0..$D){ for my $j (0..$D-$i){ next if $cur[$i][$j]->is_zero;
      for my $a (0..$D-$i-$j){ for my $b (0..$D-$i-$j-$a){ $nw[$i+$a][$j+$b] += $cur[$i][$j]*$Q[$a][$b]; } } } }
    if ($m>1) { my @n2; for my $i (0..$D){ for my $j (0..$D-$i){ $n2[$i][$j]=R(0); } }
      for my $i (0..$D){ for my $j (0..$D-$i){ next if $nw[$i][$j]->is_zero; $n2[$i+1][$j] += $nw[$i][$j] if $i+$j+1<=$D; $n2[$i][$j+1] -= $nw[$i][$j] if $i+$j+1<=$D; } }
      @nw=@n2; }
    # this is Q^m (u-u')^(m-1) only if we multiply (u-u') once per step after the first: yes
    @cur = @nw;
    for my $b (0..$mu-1){ $Am[$m][$b] = $cur[$mu-1-$b][$b] / R(fact($m)); }
  }
  my @out;
  for my $g (1..$GMAX) {
    my $E = 2*$g-1+$mu; my $tot = R(0);
    for my $b (0..$mu-1) {
      # r_b = exp( sum_k psi_k hb^k P_k(b) ), P_k(b) = sum_{c=-b}^{mu-b-1} c^k
      my $arg = ser::zero();
      for my $k (1..$E){ my $Pk=R(0); for my $c (-$b..$mu-$b-1){ $Pk += R($c)**$k; } $arg->[$k] = $psi[$k]*$Pk; }
      # exp of series with zero constant term
      my $ex = ser::zero(); $ex->[0]=R(1); my $term = ser::zero(); $term->[0]=R(1);
      for my $n (1..$E){ $term = ser::scal(R(1,$n), ser::mul($term,$arg)); $ex = ser::add($ex,$term); }
      for my $m (1..$mu){ $tot += $Am[$m][$b] * $ex->[2*$g-1+$m] if 2*$g-1+$m <= $NS; }
    }
    $tot = $tot / $mu;
    my $pred = Xi($g,$mu) * R(Math::BigInt->new($mu)**$mu) / R(fact($mu));
    push @out, "g=$g: hook=".$tot->bstr." pred=".$pred->bstr;
  }
  print "mu=$mu: ", join(" | ", @out), "\n";
}
