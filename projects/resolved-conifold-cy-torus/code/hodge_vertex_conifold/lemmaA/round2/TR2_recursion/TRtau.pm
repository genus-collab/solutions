package TRtau;
# Exact Eynard-Orantin recursion on the base curve C_tau of the dictionary P2 (Hodge CohFT, framed C^3 vertex):
#   X = log x = tau log v - (tau+1) log(1+v),  y = -v,  B = dv dv'/(v-v')^2,  one branch point v = tau.
# Local coordinate z at v = tau:  X = X(tau) - z^2/(2 sigma),  sigma = tau(1+tau);  v = tau + u(z), u = z + O(z^2) rational.
# Basis e_k(p) = Res_{q->tau} B(p,q) z(q)^(-2k-1) = Sum_{j=0}^{2k} (j+1) A_{k,j} dv/(v-tau)^(j+2),  A_{k,j} = [z^(2k)](u'(z) u(z)^j).
# omega_{g,n} = Sum C[k1..kn] prod e_{k_i}(p_i).  (Same conventions as hvc_tr_tau.wls.)
# Exact arithmetic: Math::BigRat.  ASCII only.
use strict; use warnings;
use Math::BigRat;
use Exporter 'import';
our @EXPORT_OK = qw(R new_curve run_tr xi_tr xcoef binomr moment_v);

sub R { Math::BigRat->new(@_) }
my $ZERO = R(0);

# ---------- truncated power series (non-negative powers): array refs, index = power ----------
sub ps_mul { my ($a,$b,$N)=@_; my @c=map { R(0) } 0..$N;
  for my $i (0..$#$a){ next if $i>$N; my $ai=$a->[$i]; next if $ai->is_zero;
    for my $j (0..$#$b){ last if $i+$j>$N; my $bj=$b->[$j]; next if $bj->is_zero; $c[$i+$j] += $ai*$bj; } } \@c }
sub ps_inv { my ($a,$N)=@_; my $a0=$a->[0]; die "ps_inv: zero constant" if $a0->is_zero; my @b=(R(1)/$a0);
  for my $n (1..$N){ my $s=R(0); for my $k (1..$n){ last if $k>$#$a; next if $a->[$k]->is_zero; $s += $a->[$k]*$b[$n-$k]; } push @b, -$s/$a0; } \@b }
sub ps_sqrt1 { my ($a,$N)=@_; # a[0] = 1
  my @b=(R(1)); for my $n (1..$N){ my $s=R(0); for my $k (1..$n-1){ $s += $b[$k]*$b[$n-$k]; } my $an = $n<=$#$a ? $a->[$n] : R(0); push @b, ($an-$s)/2; } \@b }
sub ps_pow { my ($a,$m,$N)=@_; my $r=[R(1), map { R(0) } 1..$N]; my $p=$a; my $e=$m;
  if ($e<0){ $p=ps_inv($a,$N); $e=-$e; }
  while ($e>0){ $r=ps_mul($r,$p,$N) if $e & 1; $e >>= 1; $p=ps_mul($p,$p,$N) if $e; } $r }
sub binomr { my ($a,$i)=@_; return R(0) if $i<0; my $r=R(1); for my $t (0..$i-1){ $r *= ($a - $t); } for my $t (1..$i){ $r /= $t; } $r }

# ---------- Laurent series: { lo => int, c => [..] } coefficient of z^(lo+i) is c->[i], truncated at power hi ----------
sub ls_new { my ($lo,$c)=@_; { lo=>$lo, c=>$c } }
sub ls_coef { my ($s,$p)=@_; my $i=$p-$s->{lo}; return $ZERO if $i<0 || $i>$#{$s->{c}}; $s->{c}[$i] }
sub ls_mul { my ($a,$b,$hi)=@_; my $lo=$a->{lo}+$b->{lo}; my $L=$hi-$lo; return ls_new($lo,[]) if $L<0;
  my @c=map { R(0) } 0..$L; my $ac=$a->{c}; my $bc=$b->{c};
  for my $i (0..$#$ac){ last if $i>$L; my $ai=$ac->[$i]; next if $ai->is_zero;
    for my $j (0..$#$bc){ last if $i+$j>$L; my $bj=$bc->[$j]; next if $bj->is_zero; $c[$i+$j] += $ai*$bj; } } ls_new($lo,\@c) }
sub ls_negz { my ($a)=@_; my @c; for my $i (0..$#{$a->{c}}){ my $p=$a->{lo}+$i; push @c, ($p % 2 == 0) ? $a->{c}[$i] : -$a->{c}[$i]; } ls_new($a->{lo},\@c) }
sub ls_add { my ($a,$b)=@_; my $lo = $a->{lo} < $b->{lo} ? $a->{lo} : $b->{lo};
  my $hia=$a->{lo}+$#{$a->{c}}; my $hib=$b->{lo}+$#{$b->{c}}; my $hi = $hia < $hib ? $hia : $hib;
  my @c; for my $p ($lo..$hi){ push @c, ls_coef($a,$p)+ls_coef($b,$p); } ls_new($lo,\@c) }

# ---------- the curve ----------
# new_curve(tau, kmax, extra): kmax = largest basis index needed (3g-3+n over the run)
sub new_curve { my ($tau,$kmax,%opt)=@_; $tau = R($tau) unless ref $tau;
  my $sig = $tau*($tau+1);
  my $HI = 2*$kmax + 4;                  # max non-negative power needed in E_k
  my $N  = 4*$kmax + 14;                 # order of u(z)
  # X(tau+u)-X(tau) = Sum_{k>=2} xk u^k,  xk = (-1)^(k+1)/k (tau^(1-k) - (tau+1)^(1-k));  z^2 = u^2 Q(u), Q_m = -2 sig x_{m+2}
  my @Q; for my $m (0..$N+2){ my $k=$m+2; my $xk = R((-1)**($k+1))/$k * ( R(1)/($tau**($k-1)) - R(1)/(($tau+1)**($k-1)) ); push @Q, -2*$sig*$xk; }
  die "Q0 != 1" unless $Q[0] == 1;
  my $Rs = ps_sqrt1(\@Q,$N+2);           # z = u R(u)
  # Lagrange: [z^n] u = (1/n) [u^(n-1)] R^(-n)
  my $Rinv = ps_inv($Rs,$N+2); my @u=(R(0)); my $P=[R(1), map { R(0) } 1..$N+2];
  for my $n (1..$N+1){ $P = ps_mul($P,$Rinv,$N+2); push @u, ($n-1 <= $#$P ? $P->[$n-1] : R(0))/$n; }
  # u(z) = Sum u[n] z^n, n=1..N+1
  my @du = map { ($_+1)*$u[$_+1] } 0..$N;          # u'(z)
  my @w = map { $u[$_+1] } 0..$N;                  # w = u/z
  my $C = { tau=>$tau, sig=>$sig, kmax=>$kmax, HI=>$HI, N=>$N, u=>\@u, du=>\@du, w=>\@w };
  # A_{k,j} = [z^(2k)](u' u^j)
  my @A; my $uj=[R(1), map { R(0) } 1..$N];  # u^0
  my @upow=($uj); my $useries=[@u]; for my $j (1..2*$kmax){ $upow[$j]=ps_mul($upow[$j-1],$useries,2*$kmax+2); }
  for my $j (0..2*$kmax){ my $p = ps_mul(\@du,$upow[$j],2*$kmax+1); for my $k (0..$kmax){ next if $j>2*$k; $A[$k][$j] = $p->[2*$k] // R(0); } }
  $C->{A}=\@A;
  # w^(-m), m = 2..2kmax+2, up to order HI + 2kmax + 2
  my $Nw = $HI + 2*$kmax + 4; my $winv = ps_inv(\@w,$Nw); my @wneg; $wneg[0]=[R(1), map { R(0) } 1..$Nw];
  for my $m (1..2*$kmax+2){ $wneg[$m]=ps_mul($wneg[$m-1],$winv,$Nw); }
  # E_k(z) = u'(z) Sum_j (j+1) A_{k,j} z^(-j-2) w^(-j-2): Laurent from z^(-2k-2) to z^HI
  my @E;
  for my $k (0..$kmax){ my $lo=-2*$k-2; my $L=$HI-$lo; my @c=map { R(0) } 0..$L;
    for my $j (0..2*$k){ my $a=$A[$k][$j]; next if $a->is_zero; my $f=$a*($j+1); my $ws=$wneg[$j+2]; my $sh=-$j-2-$lo; # z^(-j-2) offset
      for my $i (0..$L-$sh){ last if $i>$#$ws; $c[$sh+$i] += $f*$ws->[$i]; } }
    my $Ek = ls_mul(ls_new(0,[@du[0..$L]]), ls_new($lo,\@c), $HI); $E[$k]=$Ek; }
  $C->{E}=\@E;
  # Yhat = (u(z)-u(-z))/(z sig) ; invY
  my $NY = 4*$kmax+6; my @Yh; for my $i (0..$NY){ my $p=$i+1; push @Yh, ($p % 2 == 1 && $p <= $#u) ? 2*$u[$p]/$sig : R(0); }
  my $invY = ps_inv(\@Yh,$NY); $C->{invY}=$invY;
  # residue tables
  my %rhoAB; my %rhoC;
  my $iY = ls_new(0,$invY);
  my @F = map { ls_mul($E[$_],$iY,2*$kmax+2) } 0..$kmax;   # F_a = E_a invY, powers up to 2kmax+2 are needed
  my @En = map { ls_negz($E[$_]) } 0..$kmax;
  for my $a (0..$kmax){ for my $b (0..$kmax){
      my $P = ls_mul($F[$a],$En[$b],0);
      for my $k0 (0..$kmax){ my $cf = ls_coef($P,-2*$k0); next if $cf->is_zero; $rhoAB{"$k0,$a,$b"} = -$cf/(2*$k0+1); } } }
  for my $b (0..$kmax){ my $Gb = ls_mul(ls_add($E[$b],$En[$b]),$iY,0);
    for my $l (0..$b+1){ for my $k0 (0..$b+1-$l){ next if $k0>$kmax; my $cf=ls_coef($Gb,-2*$k0-2*$l); next if $cf->is_zero; $rhoC{"$k0,$l,$b"} = -$cf/(2*$k0+1); } } }
  $C->{rhoAB}=\%rhoAB; $C->{rhoC}=\%rhoC;
  # initial data
  my %W;
  $W{"0,3"} = { "0,0,0" => -2*$invY->[0] };
  # B11 = -u'(z)u'(-z)/(u(z)-u(-z))^2
  my @dun = map { ($_ % 2 == 0) ? $du[$_] : -$du[$_] } 0..$#du;
  my $num = ps_mul(\@du,\@dun,12); my @dd; for my $i (0..14){ my $p=$i+1; push @dd, ($p % 2 == 1) ? 2*$u[$p] : R(0); } # (u(z)-u(-z))/z
  my $den = ps_inv(ps_mul(\@dd,\@dd,12),12); my $b11 = ps_mul($num,$den,12); # B11 = -z^-2 b11
  my $B11 = ls_new(-2,[ map { -$_ } @$b11 ]);
  my $P11 = ls_mul($B11,$iY,0); my %w11;
  for my $k0 (0..1){ my $cf=ls_coef($P11,-2*$k0); $w11{"$k0"} = $cf/(2*$k0+1) unless $cf->is_zero; }
  $W{"1,1"}=\%w11;
  $C->{W}=\%W;
  return $C;
}

sub stable { my ($g,$n)=@_; 2*$g-2+$n > 0 }
sub addto { my ($h,$k,$v)=@_; return if $v->is_zero; if (exists $h->{$k}) { $h->{$k} += $v; delete $h->{$k} if $h->{$k}->is_zero; } else { $h->{$k} = $v->copy; } }

sub compute { my ($C,$g,$np1)=@_; my $n=$np1-1; my $W=$C->{W}; my $rhoAB=$C->{rhoAB}; my $rhoC=$C->{rhoC};
  my $k0max = 3*$g-3+$np1; my %D; my %acc;
  if ($g>=1 && stable($g-1,$n+2)) { my $Wa=$W->{($g-1).",".($n+2)} or die "missing ".($g-1).",".($n+2);
    for my $key (keys %$Wa){ addto(\%D,$key,$Wa->{$key}); } }
  for my $h (0..$g){ for my $mask (0..(2**$n)-1){ my @I=grep { ($mask>>$_)&1 } 0..$n-1; my @Ic=grep { !(($mask>>$_)&1) } 0..$n-1;
      next unless stable($h,1+@I) && stable($g-$h,1+@Ic);
      my $W1=$W->{"$h,".(1+@I)} or die "missing $h,".(1+@I); my $W2=$W->{($g-$h).",".(1+@Ic)} or die "missing ".($g-$h).",".(1+@Ic);
      for my $k1 (keys %$W1){ my @a=split /,/,$k1; my $a0=shift @a;
        for my $k2 (keys %$W2){ my @b=split /,/,$k2; my $b0=shift @b; my @kJ; my ($ia,$ib)=(0,0);
          for my $j (0..$n-1){ if (($mask>>$j)&1) { push @kJ,$a[$ia++]; } else { push @kJ,$b[$ib++]; } }
          addto(\%D, join(",",$a0,$b0,@kJ), $W1->{$k1}*$W2->{$k2}); } } } }
  for my $key (keys %D){ my @k=split /,/,$key; my $a=shift @k; my $b=shift @k; my $cf=$D{$key};
    for my $k0 (0..$k0max){ my $r=$rhoAB->{"$k0,$a,$b"} or next; addto(\%acc, join(",",$k0,@k), $cf*$r); } }
  if ($n>=1 && stable($g,$n)) { my $Wc=$W->{"$g,$n"} or die "missing $g,$n";
    for my $q (0..$n-1){ for my $key (keys %$Wc){ my @k=split /,/,$key; my $b=shift @k;
        for my $l (0..$b+1){ for my $k0 (0..$b+1-$l){ my $r=$rhoC->{"$k0,$l,$b"} or next; my @kJ=@k; splice(@kJ,$q,0,$l);
            addto(\%acc, join(",",$k0,@kJ), $Wc->{$key}*$r); } } } } }
  $W->{"$g,$np1"}=\%acc; }

# run_tr(C, gmax, nmax): all (g,n) with 2g-2+n <= 2gmax-2+nmax needed for (gmax, nmax) and below
sub run_tr { my ($C,$G,$NM,%opt)=@_; my $chimax = 2*$G-2+$NM; my @order;
  for my $chi (1..$chimax){ for my $g (0..$G){ my $n=$chi-2*$g+2; next if $n<1; next if $n > $NM + 2*($G-$g); next if ($g==0 && $n<=3) || ($g==1 && $n==1); push @order,[$g,$n]; } }
  for my $o (@order){ my $t0=time; compute($C,@$o); print STDERR "  (".join(",",@$o).") ".scalar(keys %{$C->{W}{join(",",@$o)}})." terms ".(time-$t0)."s\n" if $opt{verbose}; } }

# x-expansion at v = oo: [x^mu] xi_k = -((2k+1)/mu) Sum_{m=0}^{mu-1} [z^(2k+1)](v^(m+1)) binom((tau+1)mu, mu-m-1)
sub xicoef { my ($C,$k,$mu)=@_; my $key="$k,$mu"; return $C->{xic}{$key} if exists $C->{xic}{$key};
  my $tau=$C->{tau};
  unless ($C->{vpow}) { my $M = 2*$C->{kmax}+3; my @v = map { $_ <= $#{$C->{u}} ? $C->{u}[$_] : R(0) } 0..$M; $v[0] = $tau; $C->{vser}=\@v; $C->{vpow}=[ [R(1), map { R(0) } 1..$M] ]; }
  my $M = $#{$C->{vser}};
  while ($#{$C->{vpow}} < $mu){ push @{$C->{vpow}}, ps_mul($C->{vpow}[-1],$C->{vser},$M); }
  my $s=R(0); for my $m (0..$mu-1){ my $c = $C->{vpow}[$m+1][2*$k+1]; next if $c->is_zero; $s += $c*binomr(($tau+1)*$mu, $mu-$m-1); }
  my $r = -$s*(2*$k+1)/$mu; $C->{xic}{$key}=$r; $r }
sub xcoef { my ($C,$g,@mu)=@_; my $Wg=$C->{W}{"$g,".scalar(@mu)} or die "no W $g"; my $tot=R(0);
  for my $key (keys %$Wg){ my @k=split /,/,$key; my $p=$Wg->{$key}->copy; for my $i (0..$#mu){ $p *= xicoef($C,$k[$i],$mu[$i]); last if $p->is_zero; } $tot += $p; } $tot }
sub xi_tr { my ($C,$g,@mu)=@_; my $v = xcoef($C,$g,@mu)*($C->{sig}**($g-1)); my $tau=$C->{tau};
  for my $m (@mu){ $v /= binomr($m*($tau+1)-1,$m-1); } $v }
# Res_{v=tau} v^N e_k(v) = Sum_j (j+1) A_{k,j} binom(N, j+1) tau^(N-j-1)
sub moment_v { my ($C,$k,$N)=@_; my $tau=$C->{tau}; my $s=R(0); for my $j (0..2*$k){ next if $j+1>$N; my $a=$C->{A}[$k][$j]; next if $a->is_zero; $s += ($j+1)*$a*binomr($N,$j+1)*($tau**($N-$j-1)); } $s }
1;
