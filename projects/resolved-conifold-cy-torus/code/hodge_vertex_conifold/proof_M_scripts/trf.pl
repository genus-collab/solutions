# Float version of one-branch-point EO recursion on X = log z - z (branch point z = 1, u = z - 1); curve y chosen by name.
use strict; use warnings;
my $curve = shift // 'C'; my $GMAX = shift // 2; my $HI = shift // 26;
my $N = $HI + 8;
# --- power series (arrays 0..N) ---
sub zero { [ (0) x ($N+1) ] }
sub smul { my ($a,$b)=@_; my $c=zero(); for my $i (0..$N){ next unless $a->[$i]; for my $j (0..$N-$i){ $c->[$i+$j] += $a->[$i]*$b->[$j]; } } $c }
sub scomp { my ($f,$g)=@_; my $r=zero(); my $p=zero(); $p->[0]=1; for my $k (0..$N){ for my $i (0..$N){ $r->[$i] += $f->[$k]*$p->[$i]; } $p=smul($p,$g); } $r }
sub ssqrt1 { my ($a)=@_; my $b=zero(); $b->[0]=1; for my $n (1..$N){ my $s=0; for my $k (1..$n-1){ $s += $b->[$k]*$b->[$n-$k]; } $b->[$n]=($a->[$n]-$s)/2; } $b }
sub srevert { my ($w)=@_; my $t=zero(); $t->[1]=1; for (1..$N+1){ my $wt=scomp($w,$t); $wt->[1]-=1; for my $i (0..$N){ $t->[$i]-=$wt->[$i]; } } $t }
my $q=zero(); for my $m (0..$N){ $q->[$m] = 2*(-1)**$m/($m+2); }
my $w = ssqrt1($q); unshift @$w, 0; pop @$w;
my $T = srevert($w); my $mw = [ map { -$_ } @$w ];
my $sa = scomp($T,$mw);
# --- Laurent series: hash power => coef, keep powers <= HI ---
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
my $u = L(1=>1); my $one = L(0=>1);
my $ums = ladd($u,$s,-1);
my $G;
if ($curve eq 'H') { $G = lscal(-1, lmul(L(0=>1,1=>1), linv(lmul($u,$ums)))); }
elsif ($curve eq 'C') { my $ops=ladd($one,$s); $G = lscal(-1, lmul(lmul(L(0=>1,1=>2,2=>1),$ops), linv(lmul($u,$ums)))); }
else { die }
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
    for my $k0 (1..$K0){ my $f=ladd(L($k0=>1),$spow[$k0],-1); my $r=lres(lmul($f,$base))/2; next if abs($r)<1e-13;
      my $kk = $key eq '' ? "$k0" : "$k0,$key"; $res{$kk}+=$r; } }
  $W{"$g,".($n+1)}=\%res; }
my @order = ("0,3","1,1","1,2","2,1","0,4","1,3","2,2","3,1");
for my $gn (@order){ my ($g,$np1)=split /,/,$gn; next if $g>$GMAX; next if $g==2 && $np1==2 && $GMAX<3; next if $gn eq "0,4" && $GMAX<3; next if $gn eq "1,3" && $GMAX<3;
  my $t0=time; step($g,$np1-1); print STDERR "done $gn (",scalar(keys %{$W{$gn}})," terms) in ",time-$t0,"s\n"; }
sub fact { my $n=shift; my $f=1; $f*=$_ for 2..$n; $f }
sub binom { my ($n,$k)=@_; return 0 if $k<0; my $r=1; for my $i (1..$k){ $r=$r*($n-$k+$i)/$i; } $r }
for my $g (1..$GMAX){ my $h=$W{"$g,1"} or next; my @line;
  for my $j (1..2*$g+3){ my $tot=0;
    for my $k (keys %$h){ my $acc=0; for my $i (0..$j-1){ $acc += $j**$i/fact($i)*binom($k+$j-1-$i,$j-1-$i); } $tot += $h->{$k}*(-1)**($k+1)*$acc; }
    $tot/=$j; push @line, sprintf("Xi(%d)=%.10g", $j, $tot*fact($j)/$j**$j); }
  print "curve $curve g=$g: ", join("  ",@line), "\n"; }
print "s: ", join(", ", map { sprintf("%.6g",$sa->[$_]) } 0..6), "\n";
print "W11: ", join(", ", map { "$_ => ".sprintf("%.8g",$W{"1,1"}{$_}*24) } sort keys %{$W{"1,1"}}), "\n";
