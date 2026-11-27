use strict; use warnings; use Math::BigRat;
# Exact check: stored rows cPP_z9.m (colour 0, L = Y_0 = -alpha/(1+3t)) satisfy BS eq:pf at T_F,
# conjugated by e^{u/z} (theta u = Y): [prod_j (D+alpha) + y (3D)(3D+z)(3D+2z)] R = 0, D = Y + z theta.
# Arithmetic: Taylor series in s = t - t0 to order NS (exact BigRat). Args: file e1 e2 t0 [control]
my ($f,$e1,$e2,$t0s,$ctl)=@ARGV; $ctl//= 'none';
my $NS=3; my $R=sub{Math::BigRat->new(@_)};
local $/; open my $fh,'<',$f or die; my $s=<$fh>; close $fh;
$s=~s/\s+//g; $s=~s/^\{//; $s=~s/\}$//; my @parts; my $d=0; my $cur='';
for my $ch (split //,$s){ if($ch eq '['||$ch eq '('){$d++} elsif($ch eq ']'||$ch eq ')'){$d--}
  if($ch eq ',' && $d==0){push @parts,$cur;$cur=''} else {$cur.=$ch} } push @parts,$cur;
my $K=$#parts;
package Ser; use overload '+'=>\&add,'-'=>\&sub_,'*'=>\&mul,'/'=>\&div,'**'=>\&pw,'neg'=>\&neg;
sub new{my($c,$a)=@_; bless [@$a],$c}
sub cst{my($v)=@_; my @a=(Math::BigRat->new($v)); push @a,Math::BigRat->new(0) for 1..$NS; Ser->new(\@a)}
sub up{my($x)=@_; ref($x) eq 'Ser'?$x:cst("$x")}
sub add{my($a,$b)=@_; $b=up($b); Ser->new([map{$a->[$_]+$b->[$_]}0..$NS])}
sub sub_{my($a,$b,$sw)=@_; $b=up($b); ($a,$b)=($b,$a) if $sw; Ser->new([map{$a->[$_]-$b->[$_]}0..$NS])}
sub neg{my($a)=@_; Ser->new([map{-$a->[$_]}0..$NS])}
sub mul{my($a,$b)=@_; $b=up($b); my @c=map{Math::BigRat->new(0)}0..$NS; for my $i(0..$NS){for my $j(0..$NS-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}} Ser->new(\@c)}
sub inv{my($a)=@_; die "zero const" if $a->[0]==0; my @c=(1/$a->[0]); for my $n(1..$NS){my $v=Math::BigRat->new(0); $v+=$a->[$_]*$c[$n-$_] for 1..$n; push @c,-$v/$a->[0]} Ser->new(\@c)}
sub div{my($a,$b,$sw)=@_; $b=up($b); ($a,$b)=($b,$a) if $sw; mul($a,inv($b))}
sub pw{my($a,$b,$sw)=@_; die "swap pow" if $sw; my $n= ref($b) eq 'Ser' ? $b->[0] : Math::BigRat->new("$b"); die "nonint pow" unless $n->is_int; $n=$n->numify;
  my $r=cst(1); my $base= $n<0? inv($a):$a; for (1..abs($n)){$r=mul($r,$base)} $r}
package main;
our $t0=$R->($t0s); our ($E1,$E2)=($R->($e1),$R->($e2)); my $E=$E1+$E2; my $al=$E/3;
my @ts=($t0,$R->(1)); push @ts,$R->(0) for 2..$NS; our $T=Ser->new(\@ts);
our $LL = Ser::cst(0)-Ser::cst($al)/(Ser::cst(1)+$T*3); our $Y=$LL; our $yv=$T*$T*$T;
our ($S3,$SQ);
my $bs=chr(92);
sub ev{ my($x)=@_; my $t=$x; $t=~s/\Q$bs\E//g; $t=~s/\Q[Epsilon][1]\E/QEA/g; $t=~s/\Q[Epsilon][2]\E/QEB/g; $t=~s/\QLLL[0][y]\E/LLV/g;
  $t=~s/\QSqrt[3]\E/QSR/g; $t=~s/\bsq\b/SQV/g; $t=~s/(\d+)/Ser::cst("$1")/g; $t=~s/\^/**/g; $t=~tr/[]/()/;
  $t=~s/QEA/Ser::cst(\$E1)/g;$t=~s/QEB/Ser::cst(\$E2)/g;$t=~s/LLV/\$LL/g;$t=~s/QSR/Ser::cst(\$S3)/g;$t=~s/SQV/Ser::cst(\$SQ)/g; my $v=eval $t; die "$@ in $t" if $@; $v}
# sqrt3 exponent per entry: evaluate with S3=1 and S3=2
my @F;
for my $k(0..$K){ ($S3,$SQ)=($R->(1),$R->(1)); my $v1=ev($parts[$k]); ($S3,$SQ)=($R->(2),$R->(1)); my $v2=ev($parts[$k]);
  my $r=$v2->[0]/$v1->[0]; my $a= $r==2?1:($r==$R->(1,2)?-1:0); # entry = S3^a * v1 ; write as sqrt3 * v1 * 3^{(a-1)/2}
  my $fac = $a==1? $R->(1) : ($a==-1? $R->(1,3) : die "exp"); $F[$k]=$v1*Ser::cst($fac); }
if($ctl eq 'pert5'){ $F[5]=$F[5]*Ser::cst($R->(1001,1000)); }
if($ctl eq 'zneg'){ for my $k(0..$K){ $F[$k]=$F[$k]*Ser::cst((-1)**$k) } }
if($ctl eq 'hw'){ # multiply row by a fake z-only factor 1+z/E: must STILL satisfy PF (z-only freedom)
  my @G; for my $k(0..$K){ $G[$k]=$F[$k]; $G[$k]=$G[$k]+$F[$k-1]*Ser::cst(1/$E) if $k>0 } @F=@G; }
if($ctl eq 'tdep'){ # multiply by 1 + z t : not z-only, must FAIL
  my @G; for my $k(0..$K){ $G[$k]=$F[$k]; $G[$k]=$G[$k]+$F[$k-1]*$T if $k>0 } @F=@G; }
sub theta{my($a)=@_; my @c=map{Math::BigRat->new(0)}0..$NS; for my $n(0..$NS-1){ $c[$n]+= ($n+1)*$a->[$n+1] } my $dd=Ser->new(\@c); $dd*$T*Ser::cst($R->(1,3))}
my $M=$K+4; my $Z=sub{ [map{Ser::cst(0)}0..$M] };
sub Dop{ my($G,$c,$sh)=@_; # returns (3? ) generic: (mult*D + c z) G ; mult passed via $sh
  my @H; for my $m(0..$M){ my $v=$G->[$m]*$Y*Ser::cst($sh); $v=$v+theta($G->[$m-1])*Ser::cst($sh) if $m>0; $v=$v+$G->[$m-1]*Ser::cst($c) if $m>0 && $c!=0; $H[$m]=$v } \@H }
sub addc{ my($G,$c)=@_; [map{ $G->[$_]*Ser::cst(1) + ($c!=0? $G->[$_]*Ser::cst($c):Ser::cst(0)) }0..$M] }
my $F0=[ (map{$F[$_]}0..$K), map{Ser::cst(0)}($K+1)..$M ];
# (D+alpha)^3
my $A=$F0; for (1..3){ my $B=Dop($A,0,1); $A=[map{ $B->[$_] + $A->[$_]*Ser::cst($al) }0..$M]; }
# y (3D)(3D+z)(3D+2z)
my $Bv=$F0; $Bv=Dop($Bv,2,3); $Bv=Dop($Bv,1,3); $Bv=Dop($Bv,0,3);
my @res=map{ ($A->[$_] + $Bv->[$_]*$yv)->[0] } 0..$M;
print "e=($e1,$e2) t0=$t0s control=$ctl\n";
for my $m(0..$K+1){ my $r=$res[$m]; print "  z^$m residual: ", ($r==0?"0":"NONZERO ($r)"), "\n"; }
