#!/usr/bin/perl
# Independent check of Lemma H2 (h2 task). Own code. Rows built DIRECTLY in the true gauge from lp2_lib.m buildRmatrix:
#   P1 = C1i (L P0 + z theta P0),  P2 = (Del/C1i^2) (L P1 + z theta P1 - d2 P1),  d2 = E (C1i^2-1)/(2 Del),
# theta = y d/dy TOTAL derivative (C1i = 1/c1 a function of y), c1 = -wc * f, f = 1 + Q*Del (env Q; Q=0 is on shell).
# Normalised rows Phat = sqD P, sqD prop. to u := alpha^2 (alpha+L) (constant factors cancel against cPP's sq, Sqrt3);
# variable s = 1+3t (colour 0); L = -alpha/s. sigma := rho * h_w^{-1}(z) * u * (P0,P1,P2), rho = -s/(t t_c),
# sigma(psi) at z = -psi, lam_j = 1/(E - j psi).
use strict; use warnings; use Math::BigRat;
sub Q { Math::BigRat->new("$_[0]") }
my ($file,$E1,$E2,$NS) = @ARGV; $NS //= 10;
my $QS = Q($ENV{Q} // 0);
my ($e1,$e2) = (Q($E1),Q($E2)); my $E=$e1+$e2; my $al=$E/3;
open my $fh,'<',$file or die; my $txt = do { local $/; <$fh> }; close $fh; $txt =~ s/\s+//g; $txt =~ s/^\{//; $txt =~ s/\}$//;
my @parts; { my ($d,$c)=(0,''); for my $ch (split //,$txt){ $d++ if $ch=~/[\[\(]/; $d-- if $ch=~/[\]\)]/; if($ch eq ','&&!$d){push @parts,$c;$c=''}else{$c.=$ch} } push @parts,$c; }
my $NZ=$#parts;
our ($LV,$R3,$SQV);
sub evalm { my ($x,$Lv,$r3,$sqv)=@_; my $t=$x; $t=~s/\\//g; $t=~s/\[Epsilon\]\[1\]/(\$e1)/g; $t=~s/\[Epsilon\]\[2\]/(\$e2)/g;
  $t=~s/LLL\[0\]\[y\]/(\$LV)/g; $t=~s/Sqrt\[3\]/(\$R3)/g; $t=~s/\bsq\b/(\$SQV)/g; $t=~s/\^/**/g; $t=~tr/[]/()/;
  $t=~s/(?<![\$\w.])(\d+)/Q('$1')/g; local ($LV,$R3,$SQV)=($Lv,$r3,$sqv); my $v=eval $t; die "$@: $t" if $@; $v }
my (@pp,@Npoly);
for my $k (0..$NZ){ my $a=evalm($parts[$k],Q(7),Q(1),Q(1)); my $b=evalm($parts[$k],Q(7),Q(3),Q(1)); my $c=evalm($parts[$k],Q(7),Q(1),Q(5));
  die "sq" unless $a==$c*5; my $p = ($b/$a==3)?1:(($b/$a==Q('1/3'))?-1:die "p"); $pp[$k]=$p;
  my @X = map { Q($_)/Q(7)+Q(1) } 0..$k+3; my @Y = map { evalm($parts[$k],$_,Q(1),Q(1))*($al+$_)**($k+1) } @X;
  my @co=(Q(0))x($k+1);
  for my $i (0..$k){ my @b=(Q(1)); my $den=Q(1); for my $j (0..$k){ next if $j==$i; my @n=(Q(0))x(@b+1); for my $d (0..$#b){$n[$d+1]+=$b[$d]; $n[$d]-=$X[$j]*$b[$d];} @b=@n; $den*=($X[$i]-$X[$j]); }
    $co[$_]+=$Y[$i]*$b[$_]/$den for 0..$k; }
  for my $i ($k+1..$k+3){ my $v=Q(0); $v=$v*$X[$i]+$co[$_] for reverse 0..$k; die "interp fail k=$k" unless $v==$Y[$i]; }
  $Npoly[$k]=\@co; }
print "cPP parsed: z^0..z^$NZ, p_k = @pp; numerators (deg<=k) interpolated and verified at 3 extra points\n";
# s-series [v, [a_v, ..., a_{NS-1}]] ; WORK precision NS (exponents < NS); final readout limited by tracked loss
sub S { my ($v,@a)=@_; [$v,[@a]] }
sub g { my ($f,$e)=@_; return Q(0) if $e<$f->[0]; $f->[1][$e-$f->[0]] // Q(0) }
sub add { my ($f,$h,$c)=@_; $c//=1; my $v = $f->[0]<$h->[0]?$f->[0]:$h->[0]; [$v,[ map { g($f,$_)+$c*g($h,$_) } $v..$NS-1 ]] }
sub mul { my ($f,$h)=@_; my $v=$f->[0]+$h->[0]; my @a; for my $e ($v..$NS-1){ my $x=Q(0); for my $i ($f->[0]..$e-$h->[0]){ my $p=g($f,$i); next if $p==0; $x+=$p*g($h,$e-$i);} push @a,$x;} [$v,\@a] }
sub sc { my ($c,$f)=@_; [$f->[0],[ map { $c*$_ } @{$f->[1]} ]] }
sub inv { my ($f)=@_; my $v=$f->[0]; my @a=@{$f->[1]}; while(@a && $a[0]==0){shift @a;$v++}
  my @b=(Q(1)/$a[0]); for my $m (1..$NS-1+$v){ my $x=Q(0); for my $j (1..$m){ $x+=($a[$j]//Q(0))*$b[$m-$j]; } push @b,-$x*$b[0]; } [-$v,\@b] }
sub ds { my ($f)=@_; [$f->[0]-1,[ map { ($f->[0]+$_)*$f->[1][$_] } 0..$#{$f->[1]} ]] }
my $one=S(0,Q(1)); my $s=S(1,Q(1)); my $is=S(-1,Q(1));
my $t = S(0,Q(-1)/3,Q(1)/3);
my $th = sub { mul(sc(Q(1)/3,add($s,$one,-1)), ds($_[0])) };
my $Del = mul($s, S(0,Q(3),Q(-3),Q(1)));
my @wa=(Q(1)); push @wa, $wa[-1]*(Q($_-1)+Q(1)/3)*(Q($_-1)+Q(2)/3)/Q($_*$_) for 1..$NS+2;
my $wc=S(0,Q(1)); my $Dp=$one; my $tc=S(0,Q(0)); my @tcc=(Q(0)); { my $cum=Q(0); for my $n (1..$NS+2){ $cum+=$wa[$n-1]; push @tcc, $cum/$n; } }
for my $n (1..$NS+1){ $Dp=mul($Dp,$Del); $wc=add($wc,sc($wa[$n],$Dp)); $tc=add($tc,sc($tcc[$n],$Dp)); }
{ my $x=add($th->($tc),$wc); my $ok=1; for my $e (0..$NS-2){ $ok=0 if g($x,$e)!=0 } print "check theta t_c + wc = 0 (C1cf = -wc) to s^",$NS-2,": ",($ok?"OK":"NO"),"\n"; }
# also check wc solves theta^2 + 3y(9theta^2+9theta+2) (tcSetup L2b), y = t^3
{ my $y=mul($t,mul($t,$t)); my $a=$th->($th->($wc)); my $x=add($a, mul(sc(Q(3),$y), add(add(sc(Q(9),$a),sc(Q(9),$th->($wc))),sc(Q(2),$wc))));
  my $ok=1; for my $e (0..$NS-3){ $ok=0 if g($x,$e)!=0 } print "check L2b(wc)=0 to s^",$NS-3,": ",($ok?"OK":"NO"),"\n"; }
my $f = add($one, sc($QS,$Del));
my $c1 = sc(Q(-1), mul($wc,$f)); my $C1i = inv($c1);
my $L = sc(-$al,$is);
my $C1i2=mul($C1i,$C1i);
my $d2 = mul(sc($E/2,add($C1i2,$one,-1)), inv($Del));
my $pref2 = mul($Del, inv($C1i2));
my $u = sc($al**3, mul(add($s,$one,-1),$is));
my $iapl = inv(sc($al, mul(add($s,$one,-1),$is)));
my @Q0;
for my $k (0..$NZ){ my $num=S(0,Q(0)); my $Lp=$one; for my $j (0..$k){ $num=add($num,sc($Npoly[$k][$j],$Lp)); $Lp=mul($Lp,$L); }
  my $x=$num; $x=mul($x,$iapl) for 1..$k+1; my $e=(5+$pp[$k])/2; push @Q0, sc(Q(3)**$e,$x); }
my (@Q1,@Q2);
for my $m (0..$NZ){ my $x=mul($L,$Q0[$m]); $x=add($x,$th->($Q0[$m-1])) if $m; push @Q1, mul($C1i,$x); }
for my $m (0..$NZ){ my $x=add(mul($L,$Q1[$m]), mul($d2,$Q1[$m]),-1); $x=add($x,$th->($Q1[$m-1])) if $m; push @Q2, mul($pref2,$x); }
my @B=(Q(1)/6,Q(-1)/30,Q(1)/42,Q(-1)/30,Q(5)/66,Q(-691)/2730); my @w=($E,-$e1,-$e2); my @lg=(Q(0))x($NZ+1);
for my $k (1..6){ my $m=2*$k-1; last if $m>$NZ; my $S=Q(0); $S+=Q(1)/$_**$m for @w; $lg[$m]=$B[$k-1]/(2*$k*(2*$k-1))*$S; }
my $sgn = $ENV{HSIGN} // 1;   # h^-1 = exp(sgn * Sum lg z^m) ; edges_sing.pl convention: hinv' = Sum k lg_k hinv  (sgn=+1)
my @hi=(Q(1)); for my $n (1..$NZ){ my $v=Q(0); $v+=$_*$lg[$_]*$hi[$n-$_] for 1..$n; $hi[$n]=$sgn*$v/$n; }
my $rho = sc(Q(-1), mul(mul($s,inv($t)),inv($tc)));
my @sig;
for my $row ([@Q0],[@Q1],[@Q2]){ my @r; for my $k (0..$NZ){ my $x=S(0,Q(0)); $x=add($x,sc($hi[$k-$_],$row->[$_])) for 0..$k; push @r, mul($rho,mul($u,$x)); } push @sig,\@r; }
# lowest s-valuation used: report
{ my $mv=99; for my $a (0..2){ for my $m (0..$NZ){ my $f=$sig[$a][$m]; my $v=$f->[0]; $mv=$v if $v<$mv; } } print "min stored s-valuation of sigma entries: $mv\n"; }
my @T = map { g($tc,$_) } 1..$NS-1;
my $KT=$NS+2;
sub pmul { my ($a,$b)=@_; my @c=(Q(0))x($KT+1); for my $i (0..$#$a){ next unless $a->[$i]; for my $j (0..$#$b){ last if $i+$j>$KT; $c[$i+$j]+=$a->[$i]*($b->[$j]//Q(0)); } } \@c }
my @r=(Q(0),Q(1)/$T[0]);
for my $it (1..$KT+2){ my $sp=[@r]; my $pw=pmul($sp,$sp); my @new=(Q(0),Q(1)/$T[0],(Q(0))x($KT));
  for my $n (2..$NS-1){ my $c=$T[$n-1]//Q(0); for my $i (0..$KT){ $new[$i]-= $c*($pw->[$i]//Q(0))/$T[0]; } $pw=pmul($pw,$sp); }
  @r=@new[0..$KT]; }
my @rr = @r[1..$#r];
my @irr=(Q(1)/$rr[0]); for my $n (1..$KT-1){ my $v=Q(0); $v+=$rr[$_]*$irr[$n-$_] for 1..$n; push @irr,-$v*$irr[0]; }
my %P; $P{0}=[Q(1)]; for my $e (1..$NS){ $P{$e}=pmul($P{$e-1},\@rr) } for my $e (1..5){ $P{-$e}=pmul($P{-$e+1},\@irr) }
sub totc { my ($f,$k)=@_; my $v=Q(0); for my $e ($f->[0]..$k){ my $c=g($f,$e); next if $c==0; die "low val $e" if $e< -5; $v+=$c*($P{$e}[$k-$e]//Q(0)); } $v }
sub fit { my ($F,$J,$lin)=@_;
  my @G = map { (($_%2)?-1:1)*$F->[$_] } 0..$#$F;
  my @P=(Q(1)); for my $j (1..$J){ my @n=(Q(0))x(@P+1); for my $d (0..$#P){ $n[$d]+=$E*$P[$d]; $n[$d+1]-=$j*$P[$d]; } @P=@n; }
  my @N; for my $m (0..$#G){ my $v=Q(0); for my $d (0..($m<$#P?$m:$#P)){ $v+=$P[$d]*$G[$m-$d]; } push @N,$v; }
  my $deg=$J+($lin?1:0); return undef if $deg>=$#G;
  for my $i ($deg+1..$#G){ return undef if $N[$i]!=0 }
  my @rem=@N[0..$deg]; my @qq=(Q(0),Q(0));
  for (my $d=$deg;$d>=$J;$d--){ my $c=$rem[$d]/$P[$J]; $qq[$d-$J]=$c; $rem[$d-$J+$_]-=$c*$P[$_] for 0..$J; }
  my @beta; for my $j (1..$J){ my $x=$E/$j; my $v=Q(0); $v=$v*$x+$rem[$_] for reverse 0..$J-1; my $den=Q(1); for my $i (1..$J){ next if $i==$j; $den*=($E-$i*$x);} $beta[$j]=$v/$den; }
  return {J=>$J,a0=>$qq[0],a1=>$qq[1],beta=>\@beta,checks=>$#G-$deg}; }
my $KMAX = $ENV{KMAX} // ($NS-4);
print "e1=$e1 e2=$e2 E=$E  Q(off-shell c1=-wc(1+Q Del))=$QS  NS=$NS  readout t_c^-1..t_c^$KMAX\n";
for my $k (-1..$KMAX){ for my $a (0..2){
  my @F = map { totc($sig[$a][$_],$k) } 0..$NZ;
  my $d; for my $J (0..$NZ){ $d=fit(\@F,$J,1); last if $d; }
  my $o = sprintf("t_c^%-2d sigma_%d: ",$k,$a);
  if ($d){ my @sup=grep { $d->{beta}[$_]!=0 } 1..$d->{J}; my $lim=Q(0); $lim+=$d->{beta}[$_]/$_ for 1..$d->{J};
    $o.="J=$d->{J} checks=$d->{checks} a1=$d->{a1} a0=$d->{a0} leaves{@sup} Sum(beta/j)=$lim";
    $o.=" beta=(".join(',',map {$d->{beta}[$_]} 1..$d->{J}).")" if $k<=2; }
  else { my $d0; for my $J (0..$NZ){ $d0=fit(\@F,$J,0); last if $d0; } $o.= $d0 ? "[nolin hyp] J=$d0->{J} checks=$d0->{checks} a0=$d0->{a0} leaves{".join(' ',grep { $d0->{beta}[$_]!=0 } 1..$d0->{J})."}" : "untestable"; }
  if ($a==2){ my $pr = sc(-$E, mul(mul($Del,mul($wc,$wc)), mul($th->($f), inv($tc)))); $o.="  my_pred_a1=".totc($pr,$k); }
  print "$o\n"; } }
