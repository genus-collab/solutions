#!/usr/bin/perl
# Check that the closed form  S(s,z) = s * Sum_l r_l (s z)^l / prod_{j=1}^{l+1} (1 + j z/E)
# satisfies the T_F Picard-Fuchs operator (BS eq:pf, alpha_j = alpha = E/3) conjugated by e^{u/z}, theta u = L = 1/s - alpha:
#   q^3 Op = q^3 (1/s + z th)^3 - (N^3 + z N^2 + (2/9) z^2 N),  N = 1/s - alpha + z th,  q = 1 - alpha s,
#   th = (q/(3 alpha)) d/ds   (y = -1/(27 q^3), s = 3/(E+3L)).
# Usage: perl pfclosed.pl ALPHA K rl_out.txt [PERT_l PERT_eps]   -> exact check of z^0..z^(K+1)
use strict; use warnings; use Math::BigRat;
my ($al,$K,$rlf,$pl,$pe)=@ARGV; my $R=sub{Math::BigRat->new(@_)}; $al=$R->($al); my $E=3*$al;
open my $fh,'<',$rlf or die; my @r; while(<$fh>){ $r[$1]=$R->($2) if /r_(\d+) = (\S+)/ } close $fh;
$r[$pl]+= $R->($pe) if defined $pl;
# z-series (0..K+1) of Laurent polys in s: hash power->coef
my $Z=$K+1;
sub zero{ [ map { {} } 0..$Z ] }
sub addto{ my($A,$k,$n,$c)=@_; return if $c==0; $A->[$k]{$n}=($A->[$k]{$n}//$R->(0))+$c; }
# closed form
my $S=zero();
for my $l (0..$K){ # r_l z^l s^{l+1} * prod_j 1/(1+j z/E)
  my @t=map{$R->(0)}0..$Z; $t[$l]=$r[$l];
  for my $j(1..$l+1){ my $q=$R->($j)/$E; for my $n(1..$Z){ $t[$n]-= $q*$t[$n-1] } }
  for my $n(0..$K){ addto($S,$n,$l+1,$t[$n]) } }   # only z^0..z^K of S used
for my $n ($K+1..$Z){ $S->[$n]={} }
sub th { my($A)=@_; my $B=zero(); for my $k(0..$Z){ for my $n(keys %{$A->[$k]}){ my $c=$A->[$k]{$n};
   addto($B,$k,$n-1,$c*$n/(3*$al)); addto($B,$k,$n,-$c*$n/3); } } $B }
sub zsh { my($A,$m)=@_; my $B=zero(); for my $k(0..$Z-$m){ for my $n(keys %{$A->[$k]}){ addto($B,$k+$m,$n,$A->[$k]{$n}) } } $B }
sub smul { my($A,$p)=@_; my $B=zero(); for my $k(0..$Z){ for my $n(keys %{$A->[$k]}){ for my $e(keys %$p){ addto($B,$k,$n+$e,$A->[$k]{$n}*$p->{$e}) } } } $B }
sub add { my($A,$B,$c)=@_; $c//=1; my $C=zero(); for my $k(0..$Z){ for my $X([$A,1],[$B,$c]){ for my $n(keys %{$X->[0][$k]}){ addto($C,$k,$n,$X->[0][$k]{$n}*$X->[1]) } } } $C }
sub opM { my($A)=@_; add(smul($A,{-1=>$R->(1)}), zsh(th($A),1)) }             # (1/s + z th)
sub opN { my($A)=@_; add(smul($A,{-1=>$R->(1),0=>-$al}), zsh(th($A),1)) }       # N
my $M3=opM(opM(opM($S)));
my $q3={0=>$R->(1),1=>-3*$al,2=>3*$al**2,3=>-$al**3};
my $lhs=smul($M3,$q3);
my $N1=opN($S); my $N2=opN($N1); my $N3=opN($N2);
my $rhs=add(add($N3, zsh($N2,1)), zsh($N1,2), $R->(2,9));
my $res=add($lhs,$rhs,-1);
my $allok=1;
for my $k(0..$Z){ my @nz=grep{$res->[$k]{$_}!=0} keys %{$res->[$k]};
  printf "z^%d: %s\n",$k, @nz? "NONZERO (".scalar(@nz)." s-coefficients, e.g. s^$nz[0]: $res->[$k]{$nz[0]})" : "0"; $allok=0 if @nz; }
print $allok? "PF ANNIHILATES CLOSED FORM through z^$Z: EXACT\n" : "PF FAILS\n";
