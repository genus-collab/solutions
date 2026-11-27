# xy_point.pl: pointwise (mod p) test of formula (17) of ABDKS (Log TR through the prism of x-y swap, Prop. 3.1) for C_tau,
# n = 1, 2, 3, against the exact EO recursion (TRmod.pm).  Convention: x = X = tau log z - (tau+1) log(1+z), y = z (theirs),
#   omega_n^(g)/prod dz_i = -[hbar^(2g-2+n)] prod_i ( Sum_r T_i^r [u_i^r] ) E_1...E_n * Sum_{sigma n-cycle} prod_i 1/(w_i - wbar_sigma(i)),
#   T_i f = d/dz_i ( f / X'(z_i) ),  w_i - wbar_j = z_i - z_j + (u_i+u_j) hbar/2,  E_i = exp(-u_i (S(u_i hbar d) xhat - X)).
# usage: perl xy_point.pl tn td g n [seed]
use strict; use warnings; use integer; use FindBin; use lib $FindBin::Bin;
use TRmod qw(inv binomm);
my ($tn,$td,$g,$n,$seed)=@ARGV; $seed //= 1;
my $p = 2147483629; my $tau = TRmod::ratmod($tn,$td,$p); my $tau1=($tau+1)%$p;
srand($seed); my @zeta = map { 1000 + int(rand(1000000)) } 1..$n;
my $NH = 2*$g-2+$n;                    # hbar order
my $K = 3*$g+2;                        # Taylor truncation per variable (max derivative order r_i <= 3g+1)
my $B = $K+1; my $SZ = $B**$n;
sub idx { my @a=@_; my $s=0; my $m=1; for my $i (0..$#a){ $s += $a[$i]*$m; $m *= $B; } $s }
sub multi { my $k=shift; my @a; for (1..$n){ push @a, $k % $B; $k /= $B; } @a }
my @MI = map { [multi($_)] } 0..$SZ-1;
# univariate series helpers (length B)
sub umul { my ($a,$b)=@_; my @c=(0) x $B; for my $i (0..$K){ next unless $a->[$i]; for my $j (0..$K-$i){ next unless $b->[$j]; $c[$i+$j]=($c[$i+$j]+$a->[$i]*$b->[$j])%$p; } } \@c }
sub uinvlin { my ($c0,$e,$L)=@_; $L //= $K; # (c0 + d)^(-e) as series in d, up to d^L
  my @r; my $ic=inv($c0,$p); my $base=TRmod::powm($ic,$e,$p); for my $k (0..$L){ push @r, ($base*binomm((($p-$e)%$p),$k,$p))%$p*TRmod::powm($ic,$k,$p)%$p; } \@r }
# multivariate: hash "h,u1,..,un" -> array of length SZ
sub mv_mul_axis { my ($A,$ax,$s)=@_; my @c=(0) x $SZ; # multiply series A by univariate series s in variable ax
  my $st = $B**$ax; for my $k (0..$SZ-1){ next unless $A->[$k]; my $ai=$MI[$k][$ax]; for my $j (0..$K-$ai){ next unless $s->[$j]; my $kk=$k+$j*$st; $c[$kk]=($c[$kk]+$A->[$k]*$s->[$j])%$p; } } \@c }
sub mv_mul { my ($A,$Bm)=@_; my @c=(0) x $SZ; my @nzA = grep { $A->[$_] } 0..$SZ-1; my @nzB = grep { $Bm->[$_] } 0..$SZ-1;
  for my $i (@nzA){ my $mi=$MI[$i]; for my $j (@nzB){ my $mj=$MI[$j]; my $ok=1; my $kk=0; my $m=1; for my $a (0..$n-1){ my $s=$mi->[$a]+$mj->[$a]; if ($s>$K){ $ok=0; last; } $kk += $s*$m; $m*=$B; } next unless $ok; $c[$kk]=($c[$kk]+$A->[$i]*$Bm->[$j])%$p; } } \@c }
sub mv_add { my ($A,$Bm,$f)=@_; $f//=1; [ map { ($A->[$_] + $f*$Bm->[$_]) % $p } 0..$SZ-1 ] }
sub mv_deriv { my ($A,$ax)=@_; my @c=(0) x $SZ; my $st=$B**$ax; for my $k (0..$SZ-1){ next unless $A->[$k]; my $ai=$MI[$k][$ax]; next if $ai==0; $c[$k-$st]=($c[$k-$st]+$ai*$A->[$k])%$p; } \@c }
# 1/X'(z) = z(1+z)/(tau - z) at z = zeta + d
sub invXp { my $z=shift; my $num=[ ($z*(($z+1)%$p))%$p, (2*$z+1)%$p, 1, (0) x ($K-2) ]; my $den = uinvlin((($tau-$z)%$p+$p)%$p,1); # (tau - z - d)^(-1) = ((tau-z) - d)^-1
  # uinvlin(c0,1) gives (c0+d)^-1; we need (c0-d)^-1: flip odd signs
  my @dd = map { ($_ % 2) ? ($p-$den->[$_])%$p : $den->[$_] } 0..$K; umul($num,\@dd) }
my @IX = map { invXp($zeta[$_]) } 0..$n-1;
sub Top { my ($A,$ax)=@_; mv_deriv(mv_mul_axis($A,$ax,$IX[$ax]),$ax) }
# --- E_i: [hbar^(2k)] polynomial in u_i, univariate series in d_i ---
my @s; my @sp; { my $f=1; my @fac=(1); for my $i (1..2*$g+4){ push @fac, ($fac[-1]*$i)%$p; }
  for my $k (0..$g+1){ push @s, inv((TRmod::powm(4,$k,$p)*$fac[2*$k+1])%$p,$p); }
  @sp=(1); for my $m (1..$g+1){ my $acc=0; for my $k (1..$m){ $acc=($acc+$s[$k]*$sp[$m-$k])%$p; } push @sp, ($p-$acc)%$p; } }
my @aa=(0,$p-1); my @al=(inv($tau,$p), ($p-inv($tau1,$p))%$p);
my @Eall;   # Eall[i][k]{udeg} = univariate series
for my $i (0..$n-1){ my $z=$zeta[$i]; my @L;
  for my $k (1..$g+1){ my %P; for my $c (0..1){ my $ser = uinvlin((($z-$aa[$c])%$p+$p)%$p, 2*$k); my $f = 1; for my $t (1..2*$k-1){ $f=($f*$t)%$p; } $f=($f*inv($al[$c],$p))%$p;
      for my $pp (0..$k){ my $q=$k-$pp; my $cf = $f*$s[$pp]%$p*$sp[$q]%$p*TRmod::powm($al[$c],2*$q,$p)%$p; my $deg=2*$pp+1;
        my $old = $P{$deg} // [(0) x $B]; $P{$deg} = [ map { ($old->[$_] + $cf*$ser->[$_])%$p } 0..$K ]; } }
    $L[$k]=\%P; }
  my @E = ({ 0 => [1,(0) x $K] });
  for my $m (1..$g+1){ my %acc; for my $k (1..$m){ my $fk = ($k*inv($m,$p))%$p; for my $d1 (keys %{$L[$k]}){ for my $d2 (keys %{$E[$m-$k]}){ my $pr=umul($L[$k]{$d1},$E[$m-$k]{$d2}); my $old=$acc{$d1+$d2} // [(0) x $B]; $acc{$d1+$d2} = [ map { ($old->[$_] + $fk*$pr->[$_])%$p } 0..$K ]; } } } $E[$m]=\%acc; }
  $Eall[$i]=\@E; }
# --- cycle factor as hash "h,u_1..u_n" -> mv series ---
sub key { join(",",@_) }
my %CY;
if ($n==1) { # 1/(u hbar): represent as h=-1,u=-1
  my @one=(0) x $SZ; $one[0]=1; $CY{key(-1,-1)}=\@one; }
else {
  # all n-cycles sigma (as image lists): cyclic orders (0, c_1, ..., c_{n-1}) with (c_i) a permutation of 1..n-1
  my @perms=([]); for my $e (1..$n-1){ my @np; for my $pp (@perms){ for my $pos (0..@$pp){ my @q=@$pp; splice(@q,$pos,0,$e); push @np,\@q; } } @perms=@np; }
  my @cycles; for my $pp (@perms){ my @ord=(0,@$pp); my @img; for my $i (0..$n-1){ $img[$ord[$i]] = $ord[($i+1) % $n]; } push @cycles,\@img; }
  for my $sg (@cycles){ # product over i of 1/(z_i - z_s(i) + (u_i+u_s(i)) hbar/2) = Sum_m (-(u_i+u_s)/2)^m hbar^m (z_i-z_s)^(-m-1)
    my %acc=( key(0,(0) x $n) => do { my @o=(0) x $SZ; $o[0]=1; \@o } );
    for my $i (0..$n-1){ my $j=$sg->[$i]; my %new;
      for my $m (0..$NH){ # factor: (-1/2)^m (u_i+u_j)^m hbar^m (z_i - z_j)^(-m-1)
        my $Dz = (($zeta[$i]-$zeta[$j])%$p+$p)%$p; my $ser = uinvlin($Dz,$m+1,2*$K); # in (d_i - d_j): build mv series
        my @F=(0) x $SZ; for my $k (0..2*$K){ next unless $ser->[$k]; for my $a (0..$k){ my @mi=(0) x $n; $mi[$i]=$a; $mi[$j]=$k-$a; next if $a>$K || $k-$a>$K; my $c = binomm($k,$a,$p)*$ser->[$k]%$p; $c = ($p-$c)%$p if (($k-$a) % 2); $F[idx(@mi)] = ($F[idx(@mi)]+$c)%$p; } }
        my $hm = TRmod::powm(($p-inv(2,$p))%$p,$m,$p);
        for my $b (0..$m){ my $cb = binomm($m,$b,$p)*$hm%$p; # u_i^b u_j^(m-b)
          for my $ka (keys %acc){ my @kk=split /,/,$ka; $kk[0]+=$m; $kk[1+$i]+=$b; $kk[1+$j]+=$m-$b; next if $kk[0]>$NH; my $pr = mv_mul($acc{$ka},\@F); my $kn=key(@kk);
            $new{$kn} = mv_add($new{$kn} // [(0) x $SZ], $pr, $cb); } } }
      %acc=%new; }
    for my $ka (keys %acc){ $CY{$ka} = mv_add($CY{$ka} // [(0) x $SZ], $acc{$ka}); } } }
# --- total: -[hbar^NH] prod_i S_i  (E_1..E_n * CY) ---
my %TOT;  # key "u_1..u_n" -> mv series, at hbar order NH
for my $kc (keys %CY){ my @kk=split /,/,$kc; my $hc=shift @kk; # remaining hbar from E's: NH - hc must be even
  my $rem=$NH-$hc; next if $rem<0 || $rem%2; my $half=$rem/2;
  # distribute half = k_1+..+k_n
  my @dist=([]); for my $i (0..$n-1){ my @nd; for my $d (@dist){ my $used=0; $used+=$_ for @$d; for my $k (0..$half-$used){ push @nd,[@$d,$k]; } } @dist=@nd; }
  for my $d (@dist){ my $sum=0; $sum+=$_ for @$d; next unless $sum==$half;
    my @parts=([[ (0) x $n ], $CY{$kc}]);
    for my $i (0..$n-1){ my @np; my $Ek=$Eall[$i][$d->[$i]];
      for my $pt (@parts){ my ($us,$ser)=@$pt; for my $ud (keys %$Ek){ my @u2=@$us; $u2[$i]+=$ud; push @np, [\@u2, mv_mul_axis($ser,$i,$Ek->{$ud})]; } } @parts=@np; }
    for my $pt (@parts){ my ($us,$ser)=@$pt; my @uu = map { $us->[$_] + $kk[$_] } 0..$n-1; next if grep { $_ < 0 } @uu; my $ku=key(@uu); $TOT{$ku} = mv_add($TOT{$ku} // [(0) x $SZ], $ser); } } }
# apply S_i = Sum_r T_i^r [u_i^r] for each variable
my %cur=%TOT;
for my $i (0..$n-1){ my %next; my %groups;
  for my $ku (keys %cur){ my @uu=split /,/,$ku; my $r=$uu[$i]; $uu[$i]=0; push @{$groups{key(@uu)}}, [$r,$cur{$ku}]; }
  for my $gk (keys %groups){ my @lst = sort { $b->[0] <=> $a->[0] } @{$groups{$gk}}; my $rmax=$lst[0][0]; my %byr = map { $_->[0] => $_->[1] } @lst;
    my $acc=[(0) x $SZ]; for (my $r=$rmax; $r>=0; $r--){ $acc = mv_add($acc, $byr{$r} // [(0) x $SZ]); $acc = Top($acc,$i) if $r>0; }
    $next{$gk}=$acc; }
  %cur=%next; }
my $val = 0; for my $ku (keys %cur){ $val = ($val + $cur{$ku}[0])%$p; } $val = ($p-$val)%$p;
# --- TR side ---
my $C = TRmod::new_curve($tn,$td,3*$g-3+$n,$p); TRmod::run_tr($C,$g,$n);
my $W = $C->{W}{"$g,$n"}; my $A=$C->{A};
sub ek { my ($k,$z)=@_; my $s=0; my $iz = inv((($z-$tau)%$p+$p)%$p,$p); for my $j (0..2*$k){ my $a=$A->[$k][$j] or next; $s=($s+($j+1)*$a%$p*TRmod::powm($iz,$j+2,$p))%$p; } $s }
my $tr=0; for my $key (keys %$W){ my @k=split /,/,$key; my $q=$W->{$key}; for my $i (0..$n-1){ $q=$q*ek($k[$i],$zeta[$i])%$p; } $tr=($tr+$q)%$p; }
print "tau=$tn/$td g=$g n=$n points=@zeta: formula=$val TR=$tr ".($val==$tr ? "AGREE" : "DIFFER")."\n";
