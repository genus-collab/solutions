# Independent re-check (regular_verify).  usage: perl rv.pl entries.txt E1 E2 [t0 C1i XX PERTk PERTeps]
# 1. evaluates stored cPP_k exactly at L = 1/s - alpha for many rational s, interpolates S_k(s) := cPP_k * Sqrt3 * sq
#    (all entries normalised by the SAME prefactor Sqrt3^-1 sq^-1), and checks S_k in s Q[s], deg <= k+1 (form).
# 2. checks q^3 * [ (1/s + z th)^3 - q^-3 (N^3 + z N^2 + 2/9 z^2 N) ] S = 0 through z^9, built directly from BS eq:pf
#    on T_F (alpha_j = alpha), th = y d/dy = (q/(3 alpha)) d/ds, N = L + z th, L = 1/s - alpha, q = 1 - alpha s.
# 3. builds P0,P1,P2 exactly as lp2_lib.m buildRmatrix (with C1i, XX numeric) at t_i = w^i t0 in Q(w), and checks
#    Pi_ij(z) = F sum_{m+n=r} (-1)^n P^i_a[m] etainv^{ab} P^j_b[n] = delta_ij delta_r0, r = 0..9, F = (Sqrt3 sq)^-2.
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my ($file,$E1,$E2,$t0s,$C1is,$XXs,$pk,$pe)=@ARGV; $t0s//='1/2'; $C1is//='3/7'; $XXs//='-2/9';
my ($e1,$e2)=(R($E1),R($E2)); my $E=$e1+$e2; my $al=$E/3; my $K=9;
open my $f,'<',$file or die; my @ent=<$f>; chomp @ent; die unless @ent==10;
# ---------- parser / evaluator: value = [rat, powR, powsq] ----------
my ($str,$pos,%env);
sub peek { substr($str,$pos,1) }
sub vmul { my($x,$y)=@_; [$x->[0]*$y->[0],$x->[1]+$y->[1],$x->[2]+$y->[2]] }
sub vinv { my($x)=@_; die "div0" if $x->[0]->is_zero; [R(1)/$x->[0],-$x->[1],-$x->[2]] }
sub vadd { my($x,$y)=@_; return $y if $x->[0]->is_zero; return $x if $y->[0]->is_zero;
  die "mixed surds in sum" unless $x->[1]==$y->[1] && $x->[2]==$y->[2]; [$x->[0]+$y->[0],$x->[1],$x->[2]] }
sub vneg { my($x)=@_; [-$x->[0],$x->[1],$x->[2]] }
sub expr { my $v=term(); while($pos<length $str){ my $c=peek(); if($c eq '+'){$pos++; $v=vadd($v,term())} elsif($c eq '-'){$pos++; $v=vadd($v,vneg(term()))} else {last} } $v }
sub term { my $v=unary(); while($pos<length $str){ my $c=peek(); if($c eq '*'){$pos++; $v=vmul($v,unary())} elsif($c eq '/'){$pos++; $v=vmul($v,vinv(unary()))} elsif($c=~/[\(a-zA-Z0-9]/){ $v=vmul($v,unary()) } else {last} } $v }
sub unary { if(peek() eq '-'){ $pos++; return vneg(unary()) } if(peek() eq '+'){ $pos++; return unary() } power() }
sub power { my $b=atom(); if($pos<length $str && peek() eq '^'){ $pos++; my $neg=0; if(peek() eq '-'){$neg=1;$pos++} substr($str,$pos)=~/^(\d+)/ or die "exp"; $pos+=length $1; my $n=$1;
   my $r=[R(1),0,0]; $r=vmul($r,$b) for 1..$n; $r=vinv($r) if $neg; return $r } $b }
sub atom { my $c=peek(); if($c eq '('){ $pos++; my $v=expr(); die "paren" unless peek() eq ')'; $pos++; return $v }
  if(substr($str,$pos)=~/^(\d+)/){ $pos+=length $1; return [R($1),0,0] }
  if(substr($str,$pos,2) eq 'sq'){ $pos+=2; return [R(1),0,1] }
  if($c eq 'R'){ $pos++; return [R(1),1,0] }
  if($c=~/[abL]/){ $pos++; return [$env{$c},0,0] } die "bad token at $pos: ".substr($str,$pos,10) }
sub evalent { my($k,$L)=@_; %env=(a=>$e1,b=>$e2,L=>$L); $str=$ent[$k]; $pos=0; my $v=expr(); die "trail" if $pos!=length $str;
  die "sq power $v->[2]" unless $v->[2]==-1; my $r=$v->[0]; my $pr=$v->[1];
  # convert R^pr to R^-1 * 3^((pr+1)/2)
  die "odd" unless ($pr+1)%2==0; $r*= R(3)**(($pr+1)/2) if $pr>=-1; $r/= R(3)**(-($pr+1)/2) if $pr<-1; [$r,$pr] }
# ---------- Laurent polys in s: hash n => rat ----------
sub lp_add { my %r=%{$_[0]}; for my $n (keys %{$_[1]}){ $r{$n}=($r{$n}//R(0))+$_[1]{$n} } for(keys %r){delete $r{$_} if $r{$_}->is_zero} \%r }
sub lp_sc { my($p,$c)=@_; my %r; for(keys %$p){ my $v=$p->{$_}*$c; $r{$_}=$v unless $v->is_zero } \%r }
sub lp_mul { my($p,$q)=@_; my %r; for my $i (keys %$p){ for my $j (keys %$q){ $r{$i+$j}=($r{$i+$j}//R(0))+$p->{$i}*$q->{$j} } } for(keys %r){delete $r{$_} if $r{$_}->is_zero} \%r }
sub lp_der { my($p)=@_; my %r; for(keys %$p){ next if $_==0; $r{$_-1}=$p->{$_}*$_ } \%r }
my $Lp={-1=>R(1),0=>-$al}; my $qp={0=>R(1),1=>-$al}; my $sinv={-1=>R(1)}; my $q3=lp_mul($qp,lp_mul($qp,$qp));
sub th { lp_sc(lp_mul($qp,lp_der($_[0])), R(1)/(3*$al)) }
# ---------- 1. form ----------
my @S; my $okall=1; my @pr;
for my $k (0..$K){ my $D=$k+6; my @xs=map{R($_)} (1..$D+1+6); my @ys=map{ my $v=evalent($k,R(1)/$_ - $al); $pr[$k]=$v->[1]; $v->[0] } @xs;
  # Newton interpolation through first D+1 points, verify at remaining 6
  my @c=@ys[0..$D]; for my $j (1..$D){ for(my $i=$D;$i>=$j;$i--){ $c[$i]=($c[$i]-$c[$i-1])/($xs[$i]-$xs[$i-$j]) } }
  my @mon=(R(0)); for(my $i=$D;$i>=0;$i--){ my @n=(R(0),@mon); for my $m (0..$#mon){ $n[$m]-=$mon[$m]*$xs[$i] } $n[0]+=$c[$i]; @mon=@n }
  my $ver=1; for my $i ($D+1..$#xs){ my $v=R(0); for(my $m=$#mon;$m>=0;$m--){ $v=$v*$xs[$i]+$mon[$m] } $ver=0 unless $v==$ys[$i] }
  my @nz=grep{!$mon[$_]->is_zero}0..$#mon; my ($lo,$hi)=($nz[0],$nz[-1]);
  my $ok=$ver && $lo==1 && $hi==$k+1; $okall&&=$ok; my %h; $h{$_}=$mon[$_] for @nz; $S[$k]=\%h;
  print "form k=$k: stored Sqrt3 power $pr[$k]; interp verified on 6 extra pts: ".($ver?"yes":"NO")."; s-exponents $lo..$hi; ".($ok?"OK":"FAIL")."\n"; }
print "S_0 = $S[0]{1} s\n"; print "FORM ALL: ".($okall?"OK":"FAIL")."\n";
if(defined $pk){ $S[$pk]=lp_add($S[$pk],{1=>R($pe)}); print "CONTROL: perturbed S_$pk by $pe s\n" }
# ---------- 2. PF ----------
sub zth { my($ser)=@_; my @r=({}); for my $k (0..$K-1){ $r[$k+1]=th($ser->[$k]) } $r[$_]//={} for 0..$K; \@r }
sub sadd { my($x,$y,$c)=@_; [map{ lp_add($x->[$_], lp_sc($y->[$_],$c)) }0..$K] }
sub smulp { my($x,$p)=@_; [map{ lp_mul($x->[$_],$p) }0..$K] }
sub zsh { my($x,$m)=@_; [ (map{{}}1..$m), @{$x}[0..$K-$m] ] }
sub Aop { sadd(smulp($_[0],$sinv), zth($_[0]), R(1)) }
sub Nop { sadd(smulp($_[0],$Lp), zth($_[0]), R(1)) }
my $Sser=[@S];
my $A3=Aop(Aop(Aop($Sser))); my $N1=Nop($Sser); my $N2=Nop($N1); my $N3=Nop($N2);
my $res=smulp($A3,$q3); $res=sadd($res,$N3,R(-1)); $res=sadd($res,zsh($N2,1),R(-1)); $res=sadd($res,zsh($N1,2),R(-2)/9);
my $pfok=1; for my $j (0..$K){ my $n=scalar keys %{$res->[$j]}; $pfok=0 if $n; print "PF z^$j: ".($n?"NONZERO ($n terms)":"0")."\n" } print "PF ALL: ".($pfok?"EXACT ZERO":"FAIL")."\n";
# ---------- 3. orthogonality in Q(w) ----------
sub kn { [R($_[0]),R($_[1]//0)] } sub ka { [$_[0][0]+$_[1][0],$_[0][1]+$_[1][1]] } sub ks { [$_[0][0]*$_[1],$_[0][1]*$_[1]] }
sub km { my($x,$y)=@_; [$x->[0]*$y->[0]-$x->[1]*$y->[1], $x->[0]*$y->[1]+$x->[1]*$y->[0]-$x->[1]*$y->[1]] }
sub ki { my($x)=@_; my $N=$x->[0]**2-$x->[0]*$x->[1]+$x->[1]**2; die "kinv0" if $N->is_zero; [($x->[0]-$x->[1])/$N, -$x->[1]/$N] }
sub kz { $_[0][0]->is_zero && $_[0][1]->is_zero } sub kstr { "$_[0][0] + ($_[0][1]) w" }
sub kev { my($p,$s)=@_; my $r=kn(0); my $si=undef; for my $n (keys %$p){ my $pw=kn(1); if($n>=0){ $pw=km($pw,$s) for 1..$n } else { $si//=ki($s); $pw=km($pw,$si) for 1..-$n } $r=ka($r,ks($pw,$p->{$n})) } $r }
my $t0=R($t0s); my ($C1i,$XX)=(R($C1is),R($XXs)); my $Del=1+27*$t0**3;
my $Q1=$N1; my $NQ1=$N2;   # Q1 = N P0, NQ1 = N^2 P0  (P0 = S series)
my $k3=R(1)/(3*$e1*$e2); my @eta=([3*$al**3*$e1*$e2,0,0],[0,-$E/$k3,-1/$k3],[0,-1/$k3,0]);
my $F=$E**3/(3*$e1*$e2);
my $w=kn(0,1); my @tw=(kn($t0), km(kn($t0),$w), km(km(kn($t0),$w),$w)); my %P;
my @cols = ($t0==R(-1)/3)? (1,2) : (0,1,2);
for my $i (@cols){ my $s=km(ka(kn(1),ks($tw[$i],3)), ki(ks($tw[$i],3*$al)));
  print "colour $i: t = ".kstr($tw[$i]).", s = ".kstr($s)."\n";
  for my $k (0..$K){ my $p0=kev($S[$k],$s); my $q1=kev($Q1->[$k],$s); my $nq=kev($NQ1->[$k],$s); my $q1m=$k? kev($Q1->[$k-1],$s):kn(0);
    $P{$i}[0][$k]=$p0; $P{$i}[1][$k]=ks($q1,$C1i);
    my $in=ka(ka(ks($nq,$Del), ks($q1m,-$Del*$XX)), ks($q1,$E/2));
    $P{$i}[2][$k]=ka(ks($in,R(1)/$C1i), ks($q1,-$C1i*$E/2)); } }
my $ook=1;
for my $i (@cols){ for my $j (@cols){ my @bad;
  for my $r (0..$K){ my $acc=kn(0); for my $m (0..$r){ my $n=$r-$m; for my $a (0..2){ for my $b (0..2){ next if $eta[$a][$b]==0;
      my $v=km($P{$i}[$a][$m],$P{$j}[$b][$n]); $v=ks($v,$eta[$a][$b]*$F*($n%2?-1:1)); $acc=ka($acc,$v) } } }
    my $tgt= ($i==$j && $r==0)? kn(1):kn(0); my $d=ka($acc,ks($tgt,-1)); push @bad,"z^$r:".kstr($acc) unless kz($d) }
  $ook=0 if @bad; print "Pi_$i$j: ".(@bad? "DEVIATES @bad[0..($#bad<1?$#bad:1)]":"= delta through z^9")."\n" } }
print "ORTHO ALL: ".($ook?"EXACT":"FAIL")."\n";
