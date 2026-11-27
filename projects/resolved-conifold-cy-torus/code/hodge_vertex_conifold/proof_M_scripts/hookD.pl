# OS hook formula, s~ = delta_1 (y(z)=z): [p_mu]tau = (1/(mu mu! hb^mu)) sum_b (-1)^b C(mu-1,b) prod_{c=-b}^{mu-b-1} r(hb c), float, exact poly in hb
use strict; use warnings;
my $case = shift // 'dessin';
sub fact { my $n=shift; my $f=1; $f*=$_ for 2..$n; $f }
sub binom { my ($n,$k)=@_; return 0 if $k<0||$k>$n; my $r=1; for my $i (1..$k){ $r=$r*($n-$k+$i)/$i; } $r }
sub pmul { my ($a,$b)=@_; my @c=(0) x (@$a+@$b-1); for my $i (0..$#$a){ for my $j (0..$#$b){ $c[$i+$j]+=$a->[$i]*$b->[$j]; } } \@c }
for my $mu (1..8) { my @tot=(0) x (2*$mu+2);
  for my $b (0..$mu-1) { my $p=[1]; for my $c (-$b..$mu-$b-1) { my $f = $case eq 'dessin' ? pmul([1,$c],[1,$c]) : [1,$c]; $p = pmul($p,$f); }
    for my $k (0..$#$p){ $tot[$k] += (-1)**$b*binom($mu-1,$b)*$p->[$k]; } }
  # coefficient of hb^(2g-1) of hb^(-mu) * sum_k tot[k] hb^k  -> k = 2g-1+mu
  my @line; for my $g (1..3){ my $k=2*$g-1+$mu; my $v = $k<=$#tot ? $tot[$k] : 0; $v /= ($mu*fact($mu)); my $cm=binom(2*$mu-1,$mu-1); push @line, sprintf("g=%d: %.8g (Xi-normalised %.8g)",$g,$v,$v/$cm); }
  print "mu=$mu  ", join("  ",@line), "\n"; }
