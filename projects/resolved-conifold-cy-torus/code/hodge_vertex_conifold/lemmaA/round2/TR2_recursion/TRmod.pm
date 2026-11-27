package TRmod;
# Modular (mod p, p < 2^31 prime) version of TRtau.pm: exact Eynard-Orantin recursion on C_tau
#   X = tau log v - (tau+1) log(1+v), y = -v, B = dv dv'/(v-v')^2, branch point v = tau,
# local coordinate z: X = X(tau) - z^2/(2 sigma), v = tau + u(z); basis e_k = Sum_j (j+1) A_{k,j} dv/(v-tau)^(j+2).
# All coefficients are residues mod p (exact modular arithmetic); rational values are recovered by CRT + rational
# reconstruction over several primes (ratrec.pm).  Conventions identical to TRtau.pm / hvc_tr_tau.wls.
use strict; use warnings; use integer;
use Exporter 'import';
our @EXPORT_OK = qw(new_curve run_tr xi_tr xcoef xicoef inv mulm addm subm ratmod powm binomm);

sub inv { my ($a,$p)=@_; $a %= $p; $a += $p if $a<0; die "inv(0)" if $a==0; my ($t,$nt,$r,$nr)=(0,1,$p,$a);
  while ($nr){ my $q=$r/$nr; ($t,$nt)=($nt,$t-$q*$nt); ($r,$nr)=($nr,$r-$q*$nr); } $t += $p if $t<0; $t }
sub ratmod { my ($num,$den,$p)=@_; my $n = $num % $p; $n += $p if $n<0; my $d = $den % $p; $d += $p if $d<0; ($n*inv($d,$p)) % $p }
sub powm { my ($a,$e,$p)=@_; my $r=1; $a%=$p; $a+=$p if $a<0; if ($e<0){ $a=inv($a,$p); $e=-$e; } while ($e>0){ $r=($r*$a)%$p if $e&1; $e >>= 1; $a=($a*$a)%$p; } $r }
sub binomm { my ($a,$i,$p)=@_; return 0 if $i<0; my $r=1; for my $t (0..$i-1){ $r = ($r * (($a - $t) % $p + $p) % $p) % $p; } my $f=1; for my $t (1..$i){ $f=($f*$t)%$p; } ($r*inv($f,$p))%$p }

# power series mod p
sub ps_mul { my ($a,$b,$N,$p)=@_; my @c=(0) x ($N+1);
  for my $i (0..$#$a){ last if $i>$N; my $ai=$a->[$i] or next;
    for my $j (0..$#$b){ last if $i+$j>$N; my $bj=$b->[$j] or next; $c[$i+$j] = ($c[$i+$j] + $ai*$bj) % $p; } } \@c }
sub ps_inv { my ($a,$N,$p)=@_; my $i0=inv($a->[0],$p); my @b=($i0);
  for my $n (1..$N){ my $s=0; for my $k (1..$n){ last if $k>$#$a; my $ak=$a->[$k] or next; $s = ($s + $ak*$b[$n-$k]) % $p; } push @b, (($p-$s)*$i0) % $p; } \@b }
sub ps_sqrt1 { my ($a,$N,$p)=@_; my $i2=inv(2,$p); my @b=(1);
  for my $n (1..$N){ my $s=0; for my $k (1..$n-1){ $s=($s+$b[$k]*$b[$n-$k])%$p; } my $an = $n<=$#$a ? $a->[$n] : 0; push @b, ((($an-$s)%$p+$p)%$p*$i2)%$p; } \@b }

# Laurent series {lo, c}
sub ls_new { my ($lo,$c)=@_; { lo=>$lo, c=>$c } }
sub ls_coef { my ($s,$q)=@_; my $i=$q-$s->{lo}; return 0 if $i<0 || $i>$#{$s->{c}}; $s->{c}[$i] }
sub ls_mul { my ($a,$b,$hi,$p)=@_; my $lo=$a->{lo}+$b->{lo}; my $L=$hi-$lo; return ls_new($lo,[]) if $L<0;
  my @c=(0) x ($L+1); my $ac=$a->{c}; my $bc=$b->{c};
  for my $i (0..$#$ac){ last if $i>$L; my $ai=$ac->[$i] or next;
    for my $j (0..$#$bc){ last if $i+$j>$L; my $bj=$bc->[$j] or next; $c[$i+$j] = ($c[$i+$j] + $ai*$bj) % $p; } } ls_new($lo,\@c) }
sub ls_negz { my ($a,$p)=@_; my @c; for my $i (0..$#{$a->{c}}){ my $q=$a->{lo}+$i; my $v=$a->{c}[$i]; push @c, (($q % 2 == 0) || $v==0) ? $v : $p-$v; } ls_new($a->{lo},\@c) }
sub ls_add { my ($a,$b,$p)=@_; my $lo = $a->{lo} < $b->{lo} ? $a->{lo} : $b->{lo};
  my $hia=$a->{lo}+$#{$a->{c}}; my $hib=$b->{lo}+$#{$b->{c}}; my $hi = $hia < $hib ? $hia : $hib;
  my @c; for my $q ($lo..$hi){ push @c, (ls_coef($a,$q)+ls_coef($b,$q))%$p; } ls_new($lo,\@c) }

# new_curve(tnum, tden, kmax, p)
sub new_curve { my ($tn,$td,$kmax,$p)=@_;
  my $tau = ratmod($tn,$td,$p); my $tau1 = ($tau+1)%$p; my $sig = ($tau*$tau1)%$p;
  my $HI = 2*$kmax + 4; my $N = 4*$kmax + 14;
  my $it = inv($tau,$p); my $it1 = inv($tau1,$p);
  my @Q; for my $m (0..$N+2){ my $k=$m+2;
    my $xk = ( powm($it,$k-1,$p) - powm($it1,$k-1,$p) ) % $p; $xk += $p if $xk<0; $xk = ($xk*inv($k,$p))%$p; $xk = $p-$xk if ($k % 2 == 0) && $xk; # (-1)^(k+1)/k
    push @Q, ((($p - 2*$sig % $p) % $p) * $xk) % $p; }
  die "Q0 != 1 (got $Q[0])" unless $Q[0] == 1;
  my $Rs = ps_sqrt1(\@Q,$N+2,$p); my $Rinv = ps_inv($Rs,$N+2,$p); my @u=(0); my $P=[1,(0) x ($N+2)];
  for my $n (1..$N+1){ $P = ps_mul($P,$Rinv,$N+2,$p); push @u, (($n-1 <= $#$P ? $P->[$n-1] : 0)*inv($n,$p))%$p; }
  my @du = map { (($_+1)*$u[$_+1])%$p } 0..$N; my @w = map { $u[$_+1] } 0..$N;
  my $C = { p=>$p, tau=>$tau, tn=>$tn, td=>$td, sig=>$sig, kmax=>$kmax, HI=>$HI, N=>$N, u=>\@u, du=>\@du, w=>\@w };
  my @A; my @upow=([1,(0) x $N]); for my $j (1..2*$kmax){ $upow[$j]=ps_mul($upow[$j-1],\@u,2*$kmax+2,$p); }
  for my $j (0..2*$kmax){ my $q = ps_mul(\@du,$upow[$j],2*$kmax+1,$p); for my $k (0..$kmax){ next if $j>2*$k; $A[$k][$j] = $q->[2*$k] // 0; } }
  $C->{A}=\@A;
  my $Nw = $HI + 2*$kmax + 4; my $winv = ps_inv(\@w,$Nw,$p); my @wneg; $wneg[0]=[1,(0) x $Nw];
  for my $m (1..2*$kmax+2){ $wneg[$m]=ps_mul($wneg[$m-1],$winv,$Nw,$p); }
  my @E;
  for my $k (0..$kmax){ my $lo=-2*$k-2; my $L=$HI-$lo; my @c=(0) x ($L+1);
    for my $j (0..2*$k){ my $a=$A[$k][$j] or next; my $f=($a*($j+1))%$p; my $ws=$wneg[$j+2]; my $sh=-$j-2-$lo;
      for my $i (0..$L-$sh){ last if $i>$#$ws; $c[$sh+$i] = ($c[$sh+$i] + $f*$ws->[$i])%$p; } }
    $E[$k] = ls_mul(ls_new(0,[@du[0..$L]]), ls_new($lo,\@c), $HI, $p); }
  $C->{E}=\@E;
  my $NY = 4*$kmax+6; my $is=inv($sig,$p); my @Yh; for my $i (0..$NY){ my $q=$i+1; push @Yh, ($q % 2 == 1 && $q <= $#u) ? (2*$u[$q]%$p*$is)%$p : 0; }
  my $invY = ps_inv(\@Yh,$NY,$p); $C->{invY}=$invY;
  my %rhoAB; my %rhoC; my $iY = ls_new(0,$invY);
  my @F = map { ls_mul($E[$_],$iY,2*$kmax+2,$p) } 0..$kmax; my @En = map { ls_negz($E[$_],$p) } 0..$kmax;
  my @i2k = map { inv(2*$_+1,$p) } 0..$kmax+2;
  for my $a (0..$kmax){ for my $b (0..$kmax){ my $Pp = ls_mul($F[$a],$En[$b],0,$p);
      for my $k0 (0..$kmax){ my $cf = ls_coef($Pp,-2*$k0) or next; $rhoAB{"$k0,$a,$b"} = (($p-$cf)*$i2k[$k0])%$p; } } }
  for my $b (0..$kmax){ my $Gb = ls_mul(ls_add($E[$b],$En[$b],$p),$iY,0,$p);
    for my $l (0..$b+1){ for my $k0 (0..$b+1-$l){ next if $k0>$kmax; my $cf=ls_coef($Gb,-2*$k0-2*$l) or next; $rhoC{"$k0,$l,$b"} = (($p-$cf)*$i2k[$k0])%$p; } } }
  $C->{rhoAB}=\%rhoAB; $C->{rhoC}=\%rhoC;
  my %W; $W{"0,3"} = { "0,0,0" => ($p - (2*$invY->[0])%$p)%$p };
  my @dun = map { ($_ % 2 == 0) ? $du[$_] : ($p-$du[$_])%$p } 0..$#du;
  my $num = ps_mul(\@du,\@dun,12,$p); my @dd; for my $i (0..14){ my $q=$i+1; push @dd, ($q % 2 == 1) ? (2*$u[$q])%$p : 0; }
  my $den = ps_inv(ps_mul(\@dd,\@dd,12,$p),12,$p); my $b11 = ps_mul($num,$den,12,$p);
  my $B11 = ls_new(-2,[ map { ($p-$_)%$p } @$b11 ]); my $P11 = ls_mul($B11,$iY,0,$p); my %w11;
  for my $k0 (0..1){ my $cf=ls_coef($P11,-2*$k0) or next; $w11{"$k0"} = ($cf*$i2k[$k0])%$p; }
  $W{"1,1"}=\%w11; $C->{W}=\%W;
  return $C; }

sub stable { my ($g,$n)=@_; 2*$g-2+$n > 0 }

sub compute { my ($C,$g,$np1)=@_; my $n=$np1-1; my $p=$C->{p}; my $W=$C->{W}; my $rhoAB=$C->{rhoAB}; my $rhoC=$C->{rhoC};
  my $k0max = 3*$g-3+$np1; my %D; my %acc;
  if ($g>=1 && stable($g-1,$n+2)) { my $Wa=$W->{($g-1).",".($n+2)} or die "missing";
    while (my ($key,$v)=each %$Wa){ $D{$key} = (($D{$key}//0) + $v) % $p; } }
  for my $h (0..$g){ for my $mask (0..(2**$n)-1){ my @I=grep { ($mask>>$_)&1 } 0..$n-1; my @Ic=grep { !(($mask>>$_)&1) } 0..$n-1;
      next unless stable($h,1+@I) && stable($g-$h,1+@Ic);
      my $W1=$W->{"$h,".(1+@I)} or die "missing $h"; my $W2=$W->{($g-$h).",".(1+@Ic)} or die "missing";
      my @L1 = map { [ [split /,/,$_], $W1->{$_} ] } keys %$W1; my @L2 = map { [ [split /,/,$_], $W2->{$_} ] } keys %$W2;
      for my $e1 (@L1){ my ($a,$c1)=@$e1; for my $e2 (@L2){ my ($b,$c2)=@$e2; my @kJ; my ($ia,$ib)=(1,1);
          for my $j (0..$n-1){ if (($mask>>$j)&1) { push @kJ,$a->[$ia++]; } else { push @kJ,$b->[$ib++]; } }
          my $key=join(",",$a->[0],$b->[0],@kJ); $D{$key} = (($D{$key}//0) + $c1*$c2) % $p; } } } }
  while (my ($key,$cf)=each %D){ next unless $cf; my @k=split /,/,$key; my $a=shift @k; my $b=shift @k; my $rest=join(",",@k);
    for my $k0 (0..$k0max){ my $r=$rhoAB->{"$k0,$a,$b"} or next; my $kk = $rest eq '' ? "$k0" : "$k0,$rest"; $acc{$kk} = (($acc{$kk}//0) + $cf*$r) % $p; } }
  if ($n>=1 && stable($g,$n)) { my $Wc=$W->{"$g,$n"} or die "missing";
    for my $q (0..$n-1){ while (my ($key,$cf)=each %$Wc){ my @k=split /,/,$key; my $b=shift @k;
        for my $l (0..$b+1){ for my $k0 (0..$b+1-$l){ my $r=$rhoC->{"$k0,$l,$b"} or next; my @kJ=@k; splice(@kJ,$q,0,$l);
            my $kk=join(",",$k0,@kJ); $acc{$kk} = (($acc{$kk}//0) + $cf*$r) % $p; } } } } }
  for my $k (keys %acc){ delete $acc{$k} unless $acc{$k}; }
  $W->{"$g,$np1"}=\%acc; }

sub run_tr { my ($C,$G,$NM,%opt)=@_; my $chimax = 2*$G-2+$NM; my @order;
  for my $chi (1..$chimax){ for my $g (0..$G){ my $n=$chi-2*$g+2; next if $n<1; next if $n > $NM + 2*($G-$g); next if ($g==0 && $n<=3) || ($g==1 && $n==1); push @order,[$g,$n]; } }
  for my $o (@order){ my $t0=time; compute($C,@$o); print STDERR "  (".join(",",@$o).") ".scalar(keys %{$C->{W}{join(",",@$o)}})." terms ".(time-$t0)."s\n" if $opt{verbose}; } }

sub xicoef { my ($C,$k,$mu)=@_; my $key="$k,$mu"; return $C->{xic}{$key} if exists $C->{xic}{$key}; my $p=$C->{p}; my $tau=$C->{tau};
  unless ($C->{vpow}) { my $M = 2*$C->{kmax}+3; my @v = map { $_ <= $#{$C->{u}} ? $C->{u}[$_] : 0 } 0..$M; $v[0] = $tau; $C->{vser}=\@v; $C->{vpow}=[ [1,(0) x $M] ]; }
  my $M = $#{$C->{vser}}; while ($#{$C->{vpow}} < $mu){ push @{$C->{vpow}}, ps_mul($C->{vpow}[-1],$C->{vser},$M,$p); }
  my $a = (($tau+1)*$mu)%$p; my $s=0;
  for my $m (0..$mu-1){ my $c = $C->{vpow}[$m+1][2*$k+1] or next; $s = ($s + $c*binomm($a,$mu-$m-1,$p))%$p; }
  my $r = (($p-$s)%$p * ((2*$k+1)%$p) % $p * inv($mu,$p))%$p; $C->{xic}{$key}=$r; $r }
sub xcoef { my ($C,$g,@mu)=@_; my $p=$C->{p}; my $Wg=$C->{W}{"$g,".scalar(@mu)} or die "no W"; my $tot=0;
  while (my ($key,$v)=each %$Wg){ my @k=split /,/,$key; my $q=$v; for my $i (0..$#mu){ $q = ($q*xicoef($C,$k[$i],$mu[$i]))%$p; last unless $q; } $tot=($tot+$q)%$p; } $tot }
sub xi_tr { my ($C,$g,@mu)=@_; my $p=$C->{p}; my $v = (xcoef($C,$g,@mu)*powm($C->{sig},$g-1,$p))%$p;
  for my $m (@mu){ $v = ($v*inv(binomm((($m*($C->{tau}+1))-1)%$p,$m-1,$p),$p))%$p; } $v }
1;
