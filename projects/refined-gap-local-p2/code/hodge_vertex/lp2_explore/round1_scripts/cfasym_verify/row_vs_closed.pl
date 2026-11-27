use strict; use warnings; use Math::BigRat;
# Exact test at FINITE L (finite xi): code row cPP_z9.m, normalised by z^0 entry, times h_w^{-1}=Upsilon(-z/e1,-z/e2),
# versus BS closed form S(z,v)=Sum_l r_l (v z)^l / prod_{i=1}^{l+1}(1+i z/E), v fitted from z^1.
my ($f,$e1,$e2,$Lq,$rlfile)=@ARGV; local $/; open my $fh,'<',$f or die; my $s=<$fh>; close $fh;
$s=~s/\s+//g; $s=~s/^\{//; $s=~s/\}$//; my @parts; my $d=0; my $cur='';
for my $ch (split //,$s){ if($ch eq '['||$ch eq '('){$d++} elsif($ch eq ']'||$ch eq ')'){$d--}
  if($ch eq ',' && $d==0){push @parts,$cur;$cur=''} else {$cur.=$ch} } push @parts,$cur;
my $N=$#parts; my $R=sub{Math::BigRat->new(@_)};
our ($E1,$E2,$LL)=($R->($e1),$R->($e2),$R->($Lq)); our ($S3,$SQ);
my $bs=chr(92);
sub ev{ my($x)=@_; my $t=$x; $t=~s/\Q$bs\E//g; $t=~s/\Q[Epsilon][1]\E/QEA/g; $t=~s/\Q[Epsilon][2]\E/QEB/g; $t=~s/\QLLL[0][y]\E/LLV/g;
  $t=~s/\QSqrt[3]\E/QSR/g; $t=~s/\bsq\b/SQV/g; $t=~s/(\d+)/Math::BigRat->new("$1")/g; $t=~s/\^/**/g; $t=~tr/[]/()/;
  $t=~s/QEA/\$E1/g;$t=~s/QEB/\$E2/g;$t=~s/LLV/\$LL/g;$t=~s/QSR/\$S3/g;$t=~s/SQV/\$SQ/g; my $v=eval $t; die "$@ in $t" if $@; $v}
my (@c,@aS,@aQ);
for my $k(0..$N){ ($S3,$SQ)=($R->(1),$R->(1)); my $v1=ev($parts[$k]); ($S3,$SQ)=($R->(2),$R->(1)); my $v2=ev($parts[$k]);
  ($S3,$SQ)=($R->(1),$R->(2)); my $v3=ev($parts[$k]); my $r2=$v2/$v1; my $r3=$v3/$v1;
  my $a = $r2==2?1:$r2==$R->(1,2)?-1:$r2==1?0:die "sqrt3 exp $r2"; my $b=$r3==2?1:$r3==$R->(1,2)?-1:$r3==1?0:die "sq exp $r3";
  push @c,$v1; push @aS,$a; push @aQ,$b; }
print "sqrt3 exponents: @aS\nsq exponents: @aQ\n";
my @r; for my $k(0..$N){ my $da=$aS[$k]-$aS[0]; die "odd" if $da%2; $r[$k]=$c[$k]/$c[0]*($R->(3)**($da/2)); die "sq mismatch" if $aQ[$k]!=$aQ[0]; }
# Theta series: log Theta(x) = Sum_n B_{2n} x^{2n-1}/(2n(2n-1))
my @Bn=($R->(1,6),$R->(-1,30),$R->(1,42),$R->(-1,30),$R->(5,66),$R->(-691,2730),$R->(7,6));
my $E=$E1+$E2; my @lg=map{$R->(0)}0..$N;
# h_w^{-1} = Upsilon(-z/e1,-z/e2) = Theta(-z/e1)Theta(-z/e2)/Theta(-z/E)
for my $n(1..7){ my $m=2*$n-1; next if $m>$N; my $cf=$Bn[$n-1]/(2*$n*(2*$n-1));
  $lg[$m] += $cf*((-1/$E1)**$m + (-1/$E2)**$m - (-1/$E)**$m); }
my @h=($R->(1)); for my $n(1..$N){ my $v=$R->(0); for my $k(1..$n){$v+=$k*$lg[$k]*$h[$n-$k]} $h[$n]=$v/$n; }
sub mul{my($a,$b)=@_;my @c=map{$R->(0)}0..$N;for my $i(0..$N){for my $j(0..$N-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}}@c}
open my $rf,'<',$rlfile or die; my $rt=<$rf>; my @rl; while($rt=~/r_(\d+) = (\S+)/g){$rl[$1]=$R->($2)}
for my $sg(1,-1){ my @rs=map{$r[$_]*($sg**$_)}0..$N; my @p=mul(\@rs,\@h);
  # closed form with unknown v: [z^1] = r1 v - 1/E  => v
  my $v=($p[1]+1/$E)/$rl[1];
  my @Sx=map{$R->(0)}0..$N;
  for my $l(0..$N){ # (vz)^l * prod_{i=1}^{l+1} 1/(1+iz/E)
    my @t=map{$R->(0)}0..$N; $t[$l]=$rl[$l]*$v**$l;
    for my $i(1..$l+1){ my $q=$i/$E; for my $n(1..$N){ $t[$n]-= $q*$t[$n-1]; } } # divide by (1+qz)
    $Sx[$_]+=$t[$_] for 0..$N; }
  my @ok=map{ $p[$_]==$Sx[$_]?'=':'X' } 0..$N;
  print "row sign $sg: v=$v ; v*(E+3L)=",$v*($E+3*$LL)," ; match z^0..z^$N: @ok\n";
  print "  diff z^2: ",$p[2]-$Sx[2],"\n" if $p[2]!=$Sx[2]; }
