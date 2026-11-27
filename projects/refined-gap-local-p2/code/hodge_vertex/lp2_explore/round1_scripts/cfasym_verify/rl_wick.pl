use strict; use warnings; use Math::BigRat;
# Independent Wick computation of r_l: [exp((z/3) d1 d2) exp(s1^3/z) G(s)]_{s1=s2=0} with s3=1,
# G = (1 - X/27)^(-1/2), X = 18 s1 s2 - 4 s1^3 + s1^2 s2^2 - 4 s2^3  (from Disc/(-27 s3^2)).
my $L = shift // 8; my $A=$L; my $B=int(1.5*$L)+1;
my $R=sub{Math::BigRat->new(@_)};
# polynomials as hash "a,b" -> coeff
sub pmul { my ($p,$q)=@_; my %r; for my $k1(keys %$p){ my($a1,$b1)=split/,/,$k1; for my $k2(keys %$q){ my($a2,$b2)=split/,/,$k2;
  my($a,$b)=($a1+$a2,$b1+$b2); next if $a>$A||$b>$B; $r{"$a,$b"}=($r{"$a,$b"}//Math::BigRat->new(0))+$p->{$k1}*$q->{$k2}; }} return \%r; }
my %Y=("1,1"=>$R->(18,27),"3,0"=>$R->(-4,27),"2,2"=>$R->(1,27),"0,3"=>$R->(-4,27)); # Y = X/27
my %G=("0,0"=>$R->(1)); my $pw={"0,0"=>$R->(1)}; my $bin=$R->(1);
for my $n(1..($A+$B)){ $pw=pmul($pw,\%Y); last unless %$pw; $bin=$bin*($R->(-1,2)-$n+1)/$n; # binom(-1/2,n)
  my $c=$bin*((-1)**$n); for my $k(keys %$pw){ $G{$k}=($G{$k}//$R->(0))+$c*$pw->{$k}; } }
sub fact{my$n=shift;my$f=Math::BigRat->new(1);$f*=$_ for 1..$n;$f}
my @r=map{$R->(0)}0..$L;
for my $k(keys %G){ my($a,$b)=split/,/,$k; for(my $m=0;3*$m+$a<=$b;$m++){ next unless 3*$m+$a==$b; my $l=$b-$m; next if $l>$L;
  $r[$l]+=$G{$k}/fact($m)*fact($b)/($R->(3)**$b); } }
print "r_$_ = $r[$_]\n" for 0..$L;
