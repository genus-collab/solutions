#!/usr/bin/perl
# r_l from BS's (s1,s2)-block Wick sum (main.tex ~5257-5262), s3 = 1:
#   [exp((z/3) d1 d2) exp(s1^3/z) G(s1,s2)]_{s1=s2=0} = Sum_l r_l z^l,   G = (1-u)^{-1/2},
#   u = (s1^2 s2^2 - 4 s1^3 + 18 s1 s2 - 4 s2^3)/27   (so g(s) = -sqrt3 i s3^-2 G at s3 = 1).
# coefficient: Sum_{m,j} z^{m-j} 3^{-m} m!/j! [s1^{m-3j} s2^m] G
use strict; use warnings; use Math::BigRat;
my $LMAX = shift // 18; my $D = int(3*$LMAX/2)+1;
my $R = sub { Math::BigRat->new(@_) };
my %u = ("2,2"=>$R->(1,27), "3,0"=>$R->(-4,27), "1,1"=>$R->(18,27), "0,3"=>$R->(-4,27));
my @G; for my $a (0..$D){ for my $b (0..$D){ $G[$a][$b]=$R->(0) } }
my @P; for my $a (0..$D){ for my $b (0..$D){ $P[$a][$b]=$R->(0) } } $P[0][0]=$R->(1); # u^k
my $ck=$R->(1);
for my $k (0..2*$D){
  for my $a (0..$D){ for my $b (0..$D){ $G[$a][$b]+=$ck*$P[$a][$b] } }
  my @Q; for my $a (0..$D){ for my $b (0..$D){ $Q[$a][$b]=$R->(0) } }
  for my $a (0..$D){ for my $b (0..$D){ next if $P[$a][$b]==0;
    for my $key (keys %u){ my ($da,$db)=split /,/,$key; next if $a+$da>$D || $b+$db>$D;
      $Q[$a+$da][$b+$db]+=$P[$a][$b]*$u{$key} } } }
  @P=@Q; $ck = $ck*(2*$k+1)/(2*$k+2);   # binom(2k,k)/4^k recursion
}
my @fact=($R->(1)); $fact[$_]=$fact[$_-1]*$_ for 1..$D;
my @r=map{$R->(0)}0..$LMAX;
for my $m (0..$D){ for my $j (0..int($m/3)){ my $l=$m-$j; next if $l>$LMAX; next if $m-3*$j>$D;
  $r[$l]+= $G[$m-3*$j][$m] * $fact[$m] / $fact[$j] / ($R->(3)**$m) } }
print "r_$_ = $r[$_]\n" for 0..$LMAX;
