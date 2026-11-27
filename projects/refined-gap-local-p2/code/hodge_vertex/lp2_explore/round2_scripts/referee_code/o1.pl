use Math::BigRat; sub R{Math::BigRat->new(@_)}
for my $pt ([2,5,"1/3"],[3,-7,"-2/11"],[1,1,"5"]){
 my($a,$b,$L)=map{R($_)}@$pt; my $E=$a+$b;
 # code: cPP1/cPP0
 my $num=$a**3-6*$a**2*$b-6*$a*$b**2+$b**3+3*($a**2-11*$a*$b+$b**2)*$L;
 my $code=$num/(12*$a*$b*$E*($E+3*$L));
 # proof: omega1 + E/(12 e1 e2), c=-(E+3L)/(9L)
 my $c=-($E+3*$L)/(9*$L);
 my $pred=R(-3)/(4*$E)+R(1)/(9*$E*$c)+$E/(12*$a*$b);
 print "@$pt diff=",$code-$pred,"  (z->-z control diff=",$code-(R(-3)/(4*$E)+R(1)/(9*$E*$c)-$E/(12*$a*$b)),")\n";}
