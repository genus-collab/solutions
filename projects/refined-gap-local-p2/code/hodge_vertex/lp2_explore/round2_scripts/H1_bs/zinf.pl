#!/usr/bin/perl
# Z-infinity check: lim_{z->oo} z * (s * Sum_l r_l (s z)^l / prod_{j<=l+1}(1+jz/E)) = Sum_l r_l E^{l+1} s^{l+1}/(l+1)!
# should equal -t_c, t_c(Del) = Int_0^Del wc(x)/(1-x) dx, wc = 2F1(1/3,2/3;1;x), Del = 1 - (1 - alpha s)^-3.
use strict; use warnings; use Math::BigRat;
my ($al,$N,$rlf)=@ARGV; my $R=sub{Math::BigRat->new(@_)}; $al=$R->($al); my $E=3*$al;
open my $fh,'<',$rlf or die; my @r; while(<$fh>){ $r[$1]=$R->($2) if /r_(\d+) = (\S+)/ } close $fh;
sub mul{my($a,$b)=@_;my @c=map{$R->(0)}0..$N;for my $i(0..$N){next if $a->[$i]==0;for my $j(0..$N-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}}\@c}
# wc/(1-x) series in x, integrate
my @w=($R->(1)); for my $n(1..$N){ $w[$n]=$w[$n-1]*($R->(1,3)+$n-1)*($R->(2,3)+$n-1)/($n*$n) }
my @g=map{$R->(0)}0..$N; for my $n(0..$N){ for my $k(0..$n){ $g[$n]+=$w[$k] } }   # wc/(1-x)
my @tc=($R->(0)); for my $n(1..$N){ $tc[$n]=$g[$n-1]/$n }
# Del(s) = 1 - (1-al s)^-3 = -Sum_{n>=1} binom(n+2,2) (al s)^n
my @D=($R->(0)); for my $n(1..$N){ $D[$n]= -$R->(($n+1)*($n+2)/2)*$al**$n }
# compose tc(D(s))
my @res=map{$R->(0)}0..$N; my @P=map{$R->(0)}0..$N; $P[0]=$R->(1);
for my $k(0..$N){ for my $i(0..$N){ $res[$i]+=$tc[$k]*$P[$i] } @P=@{mul(\@P,\@D)}; }
my $f=$R->(1); my $ok=1;
for my $l(0..$N-1){ $f*= ($l+1); my $lhs=$r[$l]*$E**($l+1)/$f; my $rhs=-$res[$l+1];
  printf "s^%d: Sum-side %s ; -t_c %s ; %s\n",$l+1,$lhs,$rhs,($lhs==$rhs?"=":"X"); $ok=0 if $lhs!=$rhs; }
print $ok? "Z-INFINITY IDENTITY EXACT through s^$N\n":"MISMATCH\n";
