# Float EO recursion on X = log w - (tau+1) log(1+w), branch point w_b = 1/tau, y given by name ('base': y=-1/w).
use strict; use warnings;
my $tau = shift // 1; my $GMAX = shift // 2; my $HI = shift // 26; my $ycase = shift // 'base';
my $N = $HI + 8; my $wb = 1/$tau;
sub zero { [ (0) x ($N+1) ] }
sub smul { my ($a,$b)=@_; my $c=zero(); for my $i (0..$N){ next unless $a->[$i]; for my $j (0..$N-$i){ $c->[$i+$j] += $a->[$i]*$b->[$j]; } } $c }
sub scomp { my ($f,$g)=@_; my $r=zero(); my $p=zero(); $p->[0]=1; for my $k (0..$N){ for my $i (0..$N){ $r->[$i] += $f->[$k]*$p->[$i]; } $p=smul($p,$g); } $r }
sub ssqrt1 { my ($a)=@_; my $b=zero(); $b->[0]=1; for my $n (1..$N){ my $s=0; for my $k (1..$n-1){ $s += $b->[$k]*$b->[$n-$k]; } $b->[$n]=($a->[$n]-$s)/2; } $b }
sub srevert { my ($w)=@_; my $t=zero(); $t->[1]=1; for (1..$N+1){ my $wt=scomp($w,$t); $wt->[1]-=1; for my $i (0..$N){ $t->[$i]-=$wt->[$i]; } } $t }
# X(wb+t)-X(wb) = sum_{k>=2} Xk t^k
my @Xk = (0,0); for my $k (2..$N+2){ push @Xk, ((-1)**($k+1)/$k)*( $wb**(-$k) - ($tau+1)*(1+$wb)**(-$k) ); }
my $q=zero(); for my $m (0..$N){ $q->[$m] = $Xk[$m+2]/$Xk[2]; }
my $w = ssqrt1($q); unshift @$w, 0; pop @$w;
my $T = srevert($w); my $mw = [ map { -$_ } @$w ]; my $sa = scomp($T,$mw);
sub L { my %h=@_; \%h }
sub lmul { my ($a,$b)=@_; my %c; for my $i (keys %$a){ for my $j (keys %$b){ my $p=$i+$j; next if $p>$HI; $c{$p}+=$a->{$i}*$b->{$j}; } } \%c }
sub ladd { my ($a,$b,$sb)=@_; $sb//=1; my %c=%$a; for my $k (keys %$b){ $c{$k}+=$sb*$b->{$k}; } \%c }
sub lscal { my ($s,$a)=@_; my %c; $c{$_}=$s*$a->{$_} for keys %$a; \%c }
sub lval { my $h=shift; my @k=sort { $a <=> $b } grep { abs($h->{$_})>0 } keys %$h; $k[0] }
sub linv { my ($a)=@_; my $v=lval($a); my $a0=$a->{$v}; my $M=$HI+abs($v)+4;
  my @r=(0) x ($M+1); for my $k (keys %$a){ my $d=$k-$v; $r[$d]=$a->{$k}/$a0 if $d>=0 && $d<=$M; }
  my @b=(0) x ($M+1); $b[0]=1; for my $n (1..$M){ my $s=0; for my $k (1..$n){ $s+=$r[$k]*$b[$n-$k]; } $b[$n]=-$s; }
  my %c; for my $n (0..$M){ my $p=$n-$v; $c{$p}=$b[$n]/$a0 if $p<=$HI; } \%c }
sub fromarr { my $a=shift; my %c; for my $i (0..$#$a){ $c{$i}=$a->[$i] if $i<=$HI; } \%c }
sub lres { $_[0]->{-1} // 0 }
my $s = fromarr($sa); my @spd = map { $sa->[$_+1]*($_+1) } 0..$N-1; my $sp = fromarr(\@spd);
my $u = L(1=>1); my $one = L(0=>1); my $ums = ladd($u,$s,-1);
# X'(t) = 1/(wb+t) - (tau+1)/(1+wb+t) as series
my %Xp; for my $k (0..$HI+2){ $Xp{$k} = (-1)**$k*( $wb**(-$k-1) - ($tau+1)*(1+$wb)**(-$k-1) ); }
my $Xpser = \%Xp;
# y(t)-y(s): base y=-1/w -> (t-s)/((wb+t)(wb+s))
my $ydiff;
if ($ycase eq 'base') { $ydiff = lmul($ums, linv(lmul(L(0=>$wb,1=>1), ladd(L(0=>$wb),$s)))); }
else { die }
my $G = linv(lmul($ydiff,$Xpser));
my $AM = 6*$GMAX+4; my (@spow,@sneg,@uneg);
$spow[0]=$one; for my $k (1..$AM+1){ $spow[$k]=lmul($spow[$k-1],$s); }
for my $a (0..$AM){ $sneg[$a]=linv($spow[$a+1]); $uneg[$a]=L(-$a-1=>1); }
my %W;
sub evalfirst { my ($h,$m,$side)=@_; my %out;
  if ($h==0 && $m==1) { for my $mm (0..$AM-1){ $out{$mm+1} = $side==0 ? L($mm=>$mm+1) : lscal($mm+1, lmul($spow[$mm],$sp)); } return \%out; }
  my $Wd=$W{"$h,".(1+$m)} or return {};
  for my $key (keys %$Wd){ my @k=split /,/,$key; my $a=shift @k; my $rest=join(",",@k);
    my $t = lscal($Wd->{$key}, $side==0 ? $uneg[$a] : lmul($sneg[$a],$sp));
    $out{$rest} = exists $out{$rest} ? ladd($out{$rest},$t) : $t; } \%out }
sub mergekeys { my ($I,$ka,$kb,$n)=@_; my @a = $ka eq '' ? () : split /,/,$ka; my @b = $kb eq '' ? () : split /,/,$kb; my @k;
  for my $j (0..$n-1){ push @k, ($I->[$j] ? shift @a : shift @b); } join(",",@k) }
sub step { my ($g,$n)=@_; my %brk;
  if ($g>=1) { if ($g==1 && $n==0) { $brk{''} = lmul($sp, linv(lmul($ums,$ums))); }
    else { my $Wd=$W{($g-1).",".($n+2)}//{}; for my $key (keys %$Wd){ my @k=split /,/,$key; my $a=shift @k; my $b=shift @k; my $rest=join(",",@k);
        my $t=lscal($Wd->{$key}, lmul($uneg[$a], lmul($sneg[$b],$sp))); $brk{$rest} = exists $brk{$rest} ? ladd($brk{$rest},$t) : $t; } } }
  for my $h (0..$g){ for my $mask (0..(2**$n)-1){ my @I=map { ($mask>>$_)&1 } 0..$n-1; my $m=0; $m+=$_ for @I; my $mc=$n-$m;
      next if ($h==0 && $m==0) || ($g-$h==0 && $mc==0);
      my $A=evalfirst($h,$m,0); my $B=evalfirst($g-$h,$mc,1);
      for my $ka (keys %$A){ for my $kb (keys %$B){ my $key=mergekeys(\@I,$ka,$kb,$n); my $t=lmul($A->{$ka},$B->{$kb});
          $brk{$key} = exists $brk{$key} ? ladd($brk{$key},$t) : $t; } } } }
  my %res; my $K0=6*$g-4+2*($n+1);
  for my $key (keys %brk){ my $base=lmul($G,$brk{$key});
    for my $k0 (1..$K0){ my $f=ladd(L($k0=>1),$spow[$k0],-1); my $r=lres(lmul($f,$base))/2; next if abs($r)<1e-14;
      my $kk = $key eq '' ? "$k0" : "$k0,$key"; $res{$kk}+=$r; } }
  $W{"$g,".($n+1)}=\%res; }
our @ORDER = ("0,3","1,1","1,2","2,1"); push @ORDER, ("0,4","1,3","2,2","3,1") if $GMAX>=3;
for my $gn (@ORDER){ my ($g,$np1)=split /,/,$gn; step($g,$np1-1); }
sub binom { my ($n,$k)=@_; return 0 if $k<0; my $r=1; for my $i (1..$k){ $r=$r*($n-$k+$i)/$i; } $r }
# [x^mu]F_{g,1} = (1/mu)[w^{mu-1}] (1+w)^{mu(tau+1)} sum_k c_k (w-wb)^{-k-1};  Xi = that / binom(mu(tau+1)-1, mu-1)
sub onepoint { my ($g,$mu)=@_; my $h=$W{"$g,1"}; my $tot=0;
  for my $k (keys %$h){ my $acc=0; for my $i (0..$mu-1){ my $m=$mu-1-$i; # (1+w)^{mu(tau+1)} coeff w^i (generalized binomial)
      my $bi=1; for my $r (1..$i){ $bi *= ($mu*($tau+1)-$r+1)/$r; }
      $acc += $bi * (-1)**($k+1) * $wb**(-$k-1-$m) * binom($k+$m,$m); }
    $tot += $h->{$k}*$acc; }
  $tot /= $mu; my $cm=1; for my $r (1..$mu-1){ $cm *= ($mu*($tau+1)-$r)/$r; } $tot/$cm }
for my $g (1..$GMAX){ print "tau=$tau g=$g: ", join("  ", map { sprintf("Xi(%d)=%.10g",$_,onepoint($g,$_)) } 1..2*$g+3), "\n"; }
# expansions of omega_{g,1}/dnu at other points: nu = 1/w -> 0 (w = oo, i.e. v = 0), and rho = w+1 -> 0 (x = oo)
for my $g (1..$GMAX){ my $h=$W{"$g,1"}; my (@a,@b);
  for my $n (0..2*$g+3){ my ($sa,$sb)=(0,0);
    for my $k (keys %$h){ my $c=$h->{$k};
      # -nu^{k-1} (1 - wb nu)^{-k-1}: coefficient of nu^n is -binom(k+n-(k-1), ...) : n = k-1+m, coeff -C(k+m,m) wb^m
      my $m = $n-($k-1); $sa += -$c*binom($k+$m,$m)*$wb**$m if $m>=0;
      # 1/(rho - (1+wb))^{k+1} = (-(1+wb))^{-k-1} (1 - rho/(1+wb))^{-k-1}
      $sb += $c*(-(1+$wb))**(-$k-1)*binom($k+$n,$n)*(1+$wb)**(-$n); }
    push @a, sprintf("%.4g",$sa); push @b, sprintf("%.4g",$sb); }
  print "tau=$tau g=$g omega/dnu at w=oo: ", join(" ",@a), "\n";
  print "tau=$tau g=$g omega/drho at w=-1: ", join(" ",@b), "\n"; }
# two-point expansions at w = 0 (v = oo): omega_{g,2} = sum c_{k1,k2} prod dw_i/(w_i - wb)^{k_i+1}; coefficient of w1^{a-1} w2^{b-1} dw1 dw2
sub cw { my ($k,$p)=@_; return (-1)**($k+1)*$wb**(-$k-1-$p)*binom($k+$p,$p); }   # [w^p] (w-wb)^{-k-1}
for my $g (1..$GMAX-1){ my $h=$W{"$g,2"} or next; my @out;
  for my $a (1..2*$g+1){ for my $b ($a..2*$g+2-$a){ my $s=0; for my $key (keys %$h){ my ($k1,$k2)=split /,/,$key; $s += $h->{$key}*cw($k1,$a-1)*cw($k2,$b-1); }
      push @out, sprintf("(%d,%d):%.3g",$a,$b,$s); } }
  print "tau=$tau omega_{$g,2} w-coefficients: ", join("  ",@out), "\n"; }
