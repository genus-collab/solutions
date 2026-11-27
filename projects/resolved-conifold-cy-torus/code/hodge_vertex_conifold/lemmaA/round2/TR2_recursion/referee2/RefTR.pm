package RefTR;
# Referee-2 (TR2 route): independent Eynard-Orantin recursion on C_tau modulo a prime p, written from scratch
# (no code shared with TRmod.pm / TRtau.pm).
#   Curve: X(v) = tau log v - (tau+1) log(1+v), omega_{0,1} = -v dX, omega_{0,2} = B = dv1 dv2/(v1-v2)^2.
#   dX = (tau - v) dv/(v(1+v)): one simple branch point v = tau.
#   Local coordinate z at v = tau:  X(tau + u) - X(tau) = -z^2/(2 sigma),  sigma = tau(1+tau);  involution z -> -z.
#   Basis xi_k (k >= 0): the global 1-form with a pole only at v = tau and principal part dz/z^(2k+2) there;
#     xi_k = (1/(2k+1)) Res_{v'=tau} B(v,v') z(v')^(-2k-1) = sum_j c_{k,j} dv/(v-tau)^(j+2).
#   omega_{g,n} = sum C[k_1..k_n] prod_i xi_{k_i}(v_i)   (hash: "k1,...,kn" -> C mod p).
#   Kernel: (1/2) int_{-z}^{z} B(z0,.) = sum_k xi_k(z0) z^(2k+1);  omega_{0,1}(z) - omega_{0,1}(-z) = (v(z)-v(-z)) z dz/sigma;
#   hence omega_{g,n+1}(z0,J) = sum_k0 xi_k0(z0) (sigma/2) [z^(-2k0)] ( r(z)/Y(z) ),  Y(z) = (u(z)-u(-z))/(2z),
#   with r(z) dz^2 the usual recursion integrand (omega_{g-1,n+2}(z,-z,J) + sum' products; B-terms included).
#   Expansion at v = oo: t = 1/v, x = t (1+t)^(-tau-1);  [x^mu] int_oo xi_k = (1/mu) [t^(mu-1)] (1+t)^(mu(tau+1)) phi_k(t),
#   where xi_k = phi_k(t) dt near t = 0.
use strict; use warnings; use integer;

sub md { my ($a,$p)=@_; my $r = $a % $p; $r += $p if $r < 0; $r }
sub inv { my ($a,$p)=@_; $a = md($a,$p); die "inv(0)" if $a == 0; my ($r0,$r1,$s0,$s1)=($p,$a,0,1);
  while ($r1) { my $q = $r0 / $r1; ($r0,$r1) = ($r1, $r0 - $q*$r1); ($s0,$s1) = ($s1, $s0 - $q*$s1); } md($s0,$p) }
sub pw { my ($a,$e,$p)=@_; $a = md($a,$p); if ($e < 0) { $a = inv($a,$p); $e = -$e; } my $r = 1;
  while ($e > 0) { $r = ($r*$a) % $p if $e & 1; $a = ($a*$a) % $p; $e >>= 1; } $r }
sub rat { my ($n,$d,$p)=@_; (md($n,$p) * inv($d,$p)) % $p }
# generalized binomial C(a, l) for a mod p (a given as residue), l >= 0 integer
sub gbin { my ($a,$l,$p)=@_; return 0 if $l < 0; my $num = 1; my $den = 1;
  for my $i (0..$l-1) { $num = ($num * md($a - $i, $p)) % $p; $den = ($den * (($i+1) % $p)) % $p; } ($num * inv($den,$p)) % $p }

# ---- truncated power series (array refs, index = power, length N+1) ----
sub smul { my ($a,$b,$N,$p)=@_; my @c = (0) x ($N+1);
  for my $i (0..$#$a) { last if $i > $N; my $x = $a->[$i] or next;
    for my $j (0..$#$b) { last if $i+$j > $N; my $y = $b->[$j] or next; $c[$i+$j] = ($c[$i+$j] + $x*$y) % $p; } } \@c }
sub sinv { my ($a,$N,$p)=@_; my $i0 = inv($a->[0],$p); my @b = ($i0);
  for my $n (1..$N) { my $s = 0; for my $k (1..$n) { last if $k > $#$a; my $x = $a->[$k] or next; $s = ($s + $x*$b[$n-$k]) % $p; }
    push @b, md(-$s * $i0, $p); } \@b }
sub ssqrt1 { my ($a,$N,$p)=@_; die "sqrt: a0 != 1" unless $a->[0] == 1; my $h = inv(2,$p); my @b = (1);
  for my $n (1..$N) { my $s = 0; for my $k (1..$n-1) { $s = ($s + $b[$k]*$b[$n-$k]) % $p; }
    my $an = $n <= $#$a ? $a->[$n] : 0; push @b, (md($an - $s,$p) * $h) % $p; } \@b }
# compose sum_{n} c[n] s^n with s = series (s[0] = 0), truncated at N (Horner)
sub scomp { my ($c,$s,$N,$p)=@_; my $r = [ (0) x ($N+1) ];
  for (my $n = $#$c; $n >= 0; $n--) { $r = smul($r,$s,$N,$p); $r->[0] = ($r->[0] + ($c->[$n] // 0)) % $p; } $r }

sub new { my ($class,$tn,$td,$K,$p)=@_;
  my $S = bless { p=>$p, K=>$K, tn=>$tn, td=>$td }, $class;
  my $tau = rat($tn,$td,$p); my $t1 = ($tau+1) % $p; my $sig = ($tau*$t1) % $p;
  die "degenerate tau" if $tau == 0 || $t1 == 0;
  @$S{qw(tau t1 sig)} = ($tau,$t1,$sig);
  my $NU = 4*$K + 8;                                      # order of u(z)
  # X(tau+u) - X(tau) = sum_{k>=2} Xk u^k,  Xk = (-1)^(k-1)/k (tau^(1-k) - (1+tau)^(1-k))
  my @q; for my $m (0..$NU) { my $k = $m+2; my $xk = md(pw($tau,1-$k,$p) - pw($t1,1-$k,$p), $p);
    $xk = ($xk * inv($k,$p)) % $p; $xk = md(-$xk,$p) if $k % 2 == 0; push @q, md(-2*$sig % $p * $xk, $p); }
  die "q0 = $q[0] != 1" unless $q[0] == 1;
  my $sq = ssqrt1(\@q,$NU,$p);                              # Z(u)/u = sqrt(Q(u))
  # reversion by Lagrange: z = u/phi(u), phi = Q^(-1/2):  [z^n] u = (1/n) [u^(n-1)] phi^n
  my $phi = sinv($sq,$NU,$p); my $u = [ (0) x ($NU+1) ]; my $ph = [ 1, (0) x $NU ];
  for my $n (1..$NU) { $ph = smul($ph,$phi,$NU,$p); $u->[$n] = (($ph->[$n-1] // 0) * inv($n,$p)) % $p; }
  # check: X(tau+u(z)) - X(tau) = -z^2/(2 sigma)
  { my @Xc = (0,0, map { my $k=$_+2; md(md(pw($tau,1-$k,$p)-pw($t1,1-$k,$p),$p)*inv($k,$p)*(($k%2==0)?-1:1),$p) } 0..$NU-2);
    my $chk = scomp(\@Xc,$u,$NU,$p); my $want = md(-inv(2*$sig % $p,$p),$p);
    for my $i (0..$NU) { my $w = $i == 2 ? $want : 0; die "reversion check failed at $i" unless $chk->[$i] == $w; } }
  my @du = map { (($_+1) * $u->[$_+1]) % $p } 0..$NU-1; push @du, 0;
  $S->{u} = $u; $S->{du} = \@du; $S->{NU} = $NU;
  # A_{k,j} = [z^(2k)] (u' u^j), c_{k,j} = (j+1) A_{k,j}/(2k+1)
  my @c; my $up = [ 1, (0) x $NU ];
  my @upow = ($up); for my $j (1..2*$K) { push @upow, smul($upow[-1],$u,$NU,$p); }
  for my $j (0..2*$K) { my $q2 = smul(\@du,$upow[$j],2*$K,$p); for my $k (0..$K) { next if $j > 2*$k; $c[$k][$j] = ((($j+1) * ($q2->[2*$k] // 0)) % $p) * inv(2*$k+1,$p) % $p; } }
  $S->{c} = \@c;
  # f_k(z) = sum_j c_{k,j} u'(z) u(z)^(-(j+2)) as Laurent series; u = z w(z)
  my $MH = 2*$K + 3; my $LW = $MH + 2*$K + 4;
  my $w = [ @$u[1..$NU] ]; my $winv = sinv($w,$LW,$p);
  my @wneg = ([1,(0) x $LW]); for my $m (1..2*$K+2) { push @wneg, smul($wneg[-1],$winv,$LW,$p); }
  my @h;
  for my $k (0..$K) { my $lo = -(2*$k+2); my @L = (0) x ($MH - $lo + 1);   # index i <-> power lo + i
    for my $j (0..2*$k) { my $ckj = $c[$k][$j] or next; my $ser = smul(\@du,$wneg[$j+2],$LW,$p);   # z^(-(j+2)) * ser
      for my $i (0..$LW) { my $pow = $i - ($j+2); next if $pow < $lo || $pow > $MH; my $v = $ser->[$i] or next; $L[$pow-$lo] = ($L[$pow-$lo] + $ckj*$v) % $p; } }
    die "xi_$k: leading coefficient $L[0] != 1" unless $L[0] == 1;
    for my $i (1..2*$k+1) { die "xi_$k: non-zero singular coefficient at z^".($lo+$i) if $L[$i]; }
    $h[$k] = [ @L[(-$lo)..$#L] ]; }                       # h[k][m] = [z^m] f_k, m = 0..MH
  $S->{h} = \@h; $S->{MH} = $MH;
  # Y(z) = sum_l u_{2l+1} z^(2l), 1/Y
  my $NY = 4*$K + 8; my @Y = map { $_ % 2 == 0 && $_+1 <= $NU ? $u->[$_+1] : 0 } 0..$NY; my $iY = sinv(\@Y,$NY,$p);
  $S->{iY} = $iY;
  # rho(k0,a,b) = (sigma/2) [z^(-2k0)] ( - f_a(z) f_b(-z) / Y(z) )
  my $hs = ($sig * inv(2,$p)) % $p; my (@Pp,@Pm);
  for my $a (0..$K) { my @Ha = @{$h[$a]}; my @Hm = map { ($_ % 2) ? md(-$Ha[$_],$p) : $Ha[$_] } 0..$#Ha;
    $Pp[$a] = smul(\@Ha,$iY,$MH,$p); $Pm[$a] = smul(\@Hm,$iY,$MH,$p); }
  my @rho;
  for my $k0 (0..$K) { for my $a (0..$K) { for my $b (0..$K) { my $s = 0;
        my $e = 2*$a + 2*$b + 4 - 2*$k0; $s = ($s + $iY->[$e]) % $p if $e >= 0;
        my $e1 = 2*$a + 2 - 2*$k0; $s = ($s + $Pm[$b][$e1]) % $p if $e1 >= 0;
        my $e2 = 2*$b + 2 - 2*$k0; $s = ($s + $Pp[$a][$e2]) % $p if $e2 >= 0;
        $s = ($s + $h[$a][0]*$h[$b][0]) % $p if $k0 == 0;
        $rho[$k0][$a][$b] = md(-$hs * $s, $p); } } }
  $S->{rho} = \@rho;
  # rhoB(k0,l,b) = (sigma/2) [z^(-2k0)] ( -2(2l+1) z^(2l) f_{b,even}(z) / Y(z) )
  my @rhoB;
  for my $k0 (0..$K) { for my $b (0..$K) { for my $l (0..$b+1) { my $s = 0;
        my $e = 2*$b + 2 - 2*$k0 - 2*$l; $s = ($s + $iY->[$e]) % $p if $e >= 0;
        $s = ($s + $h[$b][0]) % $p if $k0 == 0 && $l == 0;
        $rhoB[$k0][$l][$b] = md(-$sig * ((2*$l+1) % $p) % $p * $s, $p); } } }
  $S->{rhoB} = \@rhoB;
  # base cases
  my %W; $W{"0,3"} = { "0,0,0" => md(-$sig,$p) };
  { # omega_{1,1}: r = -u'(z)u'(-z)/(4 z^2 Y^2); C[k0] = (sigma/2)[z^(-2k0)](r/Y) = -(sigma/8)[z^(2-2k0)](u'(z)u'(-z)/Y^3)
    my @dun = map { ($_ % 2) ? md(-$du[$_],$p) : $du[$_] } 0..$#du;
    my $num = smul(\@du,\@dun,4,$p); my $iY3 = smul(smul($iY,$iY,4,$p),$iY,4,$p); my $q3 = smul($num,$iY3,4,$p);
    my $f = md(-($sig * inv(8,$p)) % $p, $p); my %w11; for my $k0 (0..1) { my $v = ($f * $q3->[2-2*$k0]) % $p; $w11{"$k0"} = $v if $v; }
    $W{"1,1"} = \%w11; }
  $S->{W} = \%W;
  return $S; }

sub stable { my ($g,$n)=@_; 2*$g - 2 + $n > 0 }

sub compute { my ($S,$g,$n)=@_; my $p = $S->{p}; my $W = $S->{W}; my $rho = $S->{rho}; my $rhoB = $S->{rhoB};
  return if exists $W->{"$g,$n"};
  my $K0 = 3*$g - 3 + $n; die "K too small" if $K0 > $S->{K};
  my $m = $n - 1; my %acc;
  # term 1: omega_{g-1,n+1}(z, -z, J)
  if ($g >= 1 && stable($g-1,$n+1)) { my $Wa = $W->{($g-1).",".($n+1)} or die "missing ".($g-1).",".($n+1);
    while (my ($key,$c) = each %$Wa) { my ($a,$b,@rest) = split /,/, $key; my $rk = join(",",@rest);
      for my $k0 (0..$K0) { my $r = $rho->[$k0][$a][$b] or next; my $kk = $rk eq '' ? "$k0" : "$k0,$rk"; $acc{$kk} = (($acc{$kk} // 0) + $c*$r) % $p; } } }
  # term 2: products of two stable correlators
  for my $hh (0..$g) { for my $mask (0..(1<<$m)-1) {
      my @I = grep { ($mask >> $_) & 1 } 0..$m-1; my @Ic = grep { !(($mask >> $_) & 1) } 0..$m-1;
      my ($g1,$n1,$g2,$n2) = ($hh, 1+@I, $g-$hh, 1+@Ic);
      next unless stable($g1,$n1) && stable($g2,$n2);
      my $W1 = $W->{"$g1,$n1"} or die "missing $g1,$n1"; my $W2 = $W->{"$g2,$n2"} or die "missing $g2,$n2";
      my @L2 = map { [ [split /,/, $_], $W2->{$_} ] } keys %$W2;
      while (my ($k1,$c1) = each %$W1) { my ($a,@kI) = split /,/, $k1;
        for my $e2 (@L2) { my ($kv,$c2) = @$e2; my ($b,@kIc) = @$kv; my $cc = ($c1*$c2) % $p;
          my @out; my ($i1,$i2) = (0,0); for my $s (0..$m-1) { if (($mask >> $s) & 1) { push @out, $kI[$i1++]; } else { push @out, $kIc[$i2++]; } }
          my $rk = join(",",@out);
          for my $k0 (0..$K0) { my $r = $rho->[$k0][$a][$b] or next; my $kk = $rk eq '' ? "$k0" : "$k0,$rk"; $acc{$kk} = (($acc{$kk} // 0) + $cc*$r) % $p; } } } } }
  # term 3: B-pairs  B(z,v_j) omega_{g,n-1}(-z, J\j) + omega_{g,n-1}(z, J\j) B(-z, v_j)
  if ($m >= 1 && stable($g,$n-1)) { my $Wb = $W->{"$g,".($n-1)} or die "missing $g,".($n-1);
    for my $j (0..$m-1) { while (my ($key,$c) = each %$Wb) { my ($b,@rest) = split /,/, $key;
        for my $l (0..$b+1) { my @o = @rest; splice(@o,$j,0,$l); my $rk = join(",",@o);
          for my $k0 (0..$K0) { my $r = $rhoB->[$k0][$l][$b] or next; my $kk = "$k0,$rk"; $acc{$kk} = (($acc{$kk} // 0) + $c*$r) % $p; } } } } }
  for my $k (keys %acc) { delete $acc{$k} unless $acc{$k}; }
  $W->{"$g,$n"} = \%acc; }

# compute all (g',n') with n' <= NM + (G - g') needed for omega_{G,NM}
sub run { my ($S,$G,$NM,%o)=@_; my @list;
  for my $g (0..$G) { for my $n (1..$NM + $G - $g) { next unless stable($g,$n); next if ($g==0 && $n==3) || ($g==1 && $n==1); push @list, [$g,$n]; } }
  @list = sort { (2*$a->[0]+$a->[1]) <=> (2*$b->[0]+$b->[1]) || $a->[0] <=> $b->[0] } @list;
  for my $gn (@list) { my $t0 = time; $S->compute(@$gn); print STDERR "  W(".join(",",@$gn).") ".scalar(keys %{$S->{W}{join(",",@$gn)}})." terms, ".(time-$t0)."s\n" if $o{verbose}; } }

# X(k,mu) = [x^mu] int_oo xi_k
sub Xk { my ($S,$k,$mu)=@_; my $key = "$k,$mu"; return $S->{XC}{$key} if exists $S->{XC}{$key};
  my $p = $S->{p}; my $tau = $S->{tau}; my $N = $mu - 1;
  # phi_k(t) = - sum_j c_{k,j} t^j (1 - tau t)^(-j-2)
  my @phi = (0) x ($N+1);
  for my $j (0..2*$k) { my $ckj = $S->{c}[$k][$j] or next; next if $j > $N;
    for my $i (0..$N-$j) { # [t^i] (1 - tau t)^(-j-2) = C(j+1+i, i) tau^i
      my $cf = gbin($j+1+$i, $i, $p) * pw($tau,$i,$p) % $p; $phi[$j+$i] = md($phi[$j+$i] - $ckj*$cf, $p); } }
  my $a = ($mu * (($tau + 1) % $p)) % $p; my $s = 0;
  for my $i (0..$N) { my $f = $phi[$i] or next; $s = ($s + $f * gbin($a, $N-$i, $p)) % $p; }
  my $r = ($s * inv($mu,$p)) % $p; $S->{XC}{$key} = $r; $r }
sub xcoef { my ($S,$g,@mu)=@_; my $p = $S->{p}; my $W = $S->{W}{"$g,".scalar(@mu)} or die "no W($g,".scalar(@mu).")"; my $t = 0;
  while (my ($key,$c) = each %$W) { my @k = split /,/, $key; my $q = $c; for my $i (0..$#mu) { $q = ($q * $S->Xk($k[$i],$mu[$i])) % $p; last unless $q; } $t = ($t + $q) % $p; } $t }
# P2 normalisation: Xi_g(mu) = sigma^(g-1) [x^mu] F_{g,n} / prod C(mu_i(tau+1)-1, mu_i-1)
sub xi { my ($S,$g,@mu)=@_; my $p = $S->{p}; my $v = ($S->xcoef($g,@mu) * pw($S->{sig},$g-1,$p)) % $p;
  for my $m (@mu) { my $b = gbin(md($m*(($S->{tau}+1)%$p) - 1,$p), $m-1, $p); return undef unless $b; $v = ($v * inv($b,$p)) % $p; } $v }
1;
