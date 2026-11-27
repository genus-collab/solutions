#!/usr/bin/perl
# Symplectic pairing of the closed-form rows at one colour, all y (s symbolic), to order z^K:
#  P0 = h(z) f/sqrt(3 e1 e2 E), f = s*Sum_l r_l (s z)^l/prod_{j<=l+1}(1+jz/E);  rows by buildRmatrix (C1i, XX drop out, hand identity):
#  Pi = eta^00 P0 P0' + eta^11 P1 P1' + eta^12 (P1 P2' + P2 P1') = [alpha^3 f f' - (Q1 T' + T Q1')]/E,  ' = (z -> -z),
#  Q1 = N f, T = Del N Q1 + (E/2) Q1, N = 1/s - alpha + z th, th = (q/(3 alpha)) d/ds, q = 1 - alpha s, Del = 1 - q^-3.
#  Checked in the form  alpha^3 q^3 f f' - (Q1 That' + That Q1') - E q^3 = 0,  That = q^3 T = (q^3-1) N Q1 + (E/2) q^3 Q1.
# Usage: perl symp.pl ALPHA K rl_out.txt [PERT_l eps]
use strict; use warnings; use Math::BigRat;
my ($al,$K,$rlf,$pl,$pe)=@ARGV; my $R=sub{Math::BigRat->new(@_)}; $al=$R->($al); my $E=3*$al; my $Z=$K;
open my $fh,'<',$rlf or die; my @r; while(<$fh>){ $r[$1]=$R->($2) if /r_(\d+) = (\S+)/ } close $fh;
$r[$pl]+=$R->($pe) if defined $pl;
sub zero{ [ map { {} } 0..$Z ] }
sub addto{ my($A,$k,$n,$c)=@_; return if $c==0; $A->[$k]{$n}=($A->[$k]{$n}//$R->(0))+$c; }
my $f=zero();
for my $l (0..$K){ my @t=map{$R->(0)}0..$Z; $t[$l]=$r[$l];
  for my $j(1..$l+1){ my $q=$R->($j)/$E; for my $n(1..$Z){ $t[$n]-= $q*$t[$n-1] } }
  for my $n(0..$Z){ addto($f,$n,$l+1,$t[$n]) } }
sub th { my($A)=@_; my $B=zero(); for my $k(0..$Z){ for my $n(keys %{$A->[$k]}){ my $c=$A->[$k]{$n};
   addto($B,$k,$n-1,$c*$n/(3*$al)); addto($B,$k,$n,-$c*$n/3); } } $B }
sub zsh { my($A,$m)=@_; my $B=zero(); for my $k(0..$Z-$m){ for my $n(keys %{$A->[$k]}){ addto($B,$k+$m,$n,$A->[$k]{$n}) } } $B }
sub smul { my($A,$p)=@_; my $B=zero(); for my $k(0..$Z){ for my $n(keys %{$A->[$k]}){ for my $e(keys %$p){ addto($B,$k,$n+$e,$A->[$k]{$n}*$p->{$e}) } } } $B }
sub add { my($A,$B,$c)=@_; $c//=1; my $C=zero(); for my $k(0..$Z){ for my $X([$A,1],[$B,$c]){ for my $n(keys %{$X->[0][$k]}){ addto($C,$k,$n,$X->[0][$k]{$n}*$X->[1]) } } } $C }
sub flip { my($A)=@_; my $B=zero(); for my $k(0..$Z){ for my $n(keys %{$A->[$k]}){ addto($B,$k,$n,$A->[$k]{$n}*((-1)**$k)) } } $B }
sub mul { my($A,$B)=@_; my $C=zero(); for my $i(0..$Z){ for my $j(0..$Z-$i){ for my $n(keys %{$A->[$i]}){ for my $m(keys %{$B->[$j]}){ addto($C,$i+$j,$n+$m,$A->[$i]{$n}*$B->[$j]{$m}) } } } } $C }
sub opN { my($A)=@_; add(smul($A,{-1=>$R->(1),0=>-$al}), zsh(th($A),1)) }
my $q3={0=>$R->(1),1=>-3*$al,2=>3*$al**2,3=>-$al**3}; my $q3m1={1=>-3*$al,2=>3*$al**2,3=>-$al**3};
my $Q1=opN($f); my $That=add(smul(opN($Q1),$q3m1), smul($Q1,$q3), $E/2);
my $lhs=add( smul(mul($f,flip($f)),{map{$_=>$q3->{$_}*$al**3} keys %$q3}),
             add(mul($Q1,flip($That)), mul($That,flip($Q1))), -1);
my $res=add($lhs, do{ my $c=zero(); for my $e(keys %$q3){ addto($c,0,$e,$E*$q3->{$e}) } $c }, -1);
my $ok=1; for my $k(0..$Z){ my @nz=grep{$res->[$k]{$_}!=0} keys %{$res->[$k]}; printf "z^%d: %s\n",$k,@nz?"NONZERO":"0"; $ok=0 if @nz }
print $ok? "Pi_ii = 1 identically in s through z^$Z: EXACT\n" : "Pi_ii != 1\n";
