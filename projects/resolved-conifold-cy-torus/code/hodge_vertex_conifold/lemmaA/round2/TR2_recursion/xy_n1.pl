# xy_n1.pl: the one-point function of C_tau from the closed formula of Alexandrov-Bychkov-Dunin-Barkowski-Kazarian-Shadrin,
# "Log topological recursion through the prism of x-y swap" (IMRN 2024, arXiv:2312.16950), Proposition 3.1 (y affine), n = 1:
#   omega_{g,1} = [hbar^{2g}] Sum_{r>=0} (d (1/dX))^r [u^{r+1}] E  dv ,   (our y = -v; the sign of w - wbar is -u hbar)
#   E = exp(-u (S(u hbar d_v) xhat - X)),  xhat = Sum_i (1/S(alpha_i hbar d_v)) log(v - a_i)/alpha_i,  S(t) = 2 sinh(t/2)/t,
#   a_0 = 0, alpha_0 = 1/tau;  a_1 = -1, alpha_1 = -1/(tau+1)   (1/alpha_i = residues of dX at the LogTR-vital points of the dual).
# Everything is expanded at v = oo in t = 1/v; output: Xi_g(mu) = sigma^(g-1) [x^mu] F_{g,1} / binom(mu(tau+1)-1, mu-1)  (P2 normalisation).
# usage: perl xy_n1.pl tau G M      (exact, Math::BigRat; tau rational, e.g. 2 or 3/7)
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $tau = R(shift // 2); my $G = shift // 3; my $M = shift // (2*$G+2);
my $sig = $tau*($tau+1);
my $T = $M + 4;                                 # t-series truncation order
sub smul { my ($a,$b)=@_; my @c=map { R(0) } 0..$T; for my $i (0..$T){ next if $a->[$i]->is_zero; for my $j (0..$T-$i){ next if $b->[$j]->is_zero; $c[$i+$j] += $a->[$i]*$b->[$j]; } } \@c }
sub sadd { my ($a,$b,$s)=@_; $s //= R(1); [ map { $a->[$_] + $s*$b->[$_] } 0..$T ] }
sub sscal { my ($s,$a)=@_; [ map { $s*$a->[$_] } 0..$T ] }
sub zero { [ map { R(0) } 0..$T ] }
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
sub binomr { my ($a,$i)=@_; return R(0) if $i<0; my $r=R(1); for my $k (0..$i-1){ $r *= ($a-$k); } $r/fact($i) }
# S(t) = Sum s_k t^(2k), s_k = 1/(4^k (2k+1)!);  1/S = Sum sp_k t^(2k)
my $K = $G + 1; my @s = map { R(1)/(R(4)**$_ * fact(2*$_+1)) } 0..$K;
my @sp = (R(1)); for my $n (1..$K){ my $acc=R(0); for my $k (1..$n){ $acc += $s[$k]*$sp[$n-$k]; } push @sp, -$acc; }
my @a = (R(0), R(-1)); my @al = (R(1)/$tau, R(-1)/($tau+1));
# (v - a)^(-m) as a t-series: t^m (1 - a t)^(-m)
sub vpow_inv { my ($aa,$m)=@_; my $r=zero(); for my $i (0..$T-$m){ $r->[$m+$i] = binomr(R(-$m),$i)*((-$aa)**$i); } $r }
# L_k(u,t) = u (2k-1)! Sum_i c_k(u,alpha_i)/(alpha_i (v-a_i)^(2k)),  c_k(u,al) = Sum_{p+q=k} s_p sp_q u^(2p) al^(2q)
# store polynomials in u as hash power -> t-series
my @L; for my $k (1..$G){ my %P;
  for my $i (0..1){ my $vi = vpow_inv($a[$i],2*$k); my $f = fact(2*$k-1)/$al[$i];
    for my $p (0..$k){ my $q=$k-$p; my $cf = $f*$s[$p]*$sp[$q]*($al[$i]**(2*$q)); my $deg = 2*$p+1;
      $P{$deg} = exists $P{$deg} ? sadd($P{$deg},$vi,$cf) : sscal($cf,$vi); } }
  $L[$k]=\%P; }
sub pmul { my ($A,$B)=@_; my %C; for my $i (keys %$A){ for my $j (keys %$B){ my $t=smul($A->{$i},$B->{$j}); $C{$i+$j} = exists $C{$i+$j} ? sadd($C{$i+$j},$t) : $t; } } \%C }
sub padd { my ($A,$B,$s)=@_; my %C = map { $_ => $A->{$_} } keys %$A; for my $j (keys %$B){ $C{$j} = exists $C{$j} ? sadd($C{$j},$B->{$j},$s) : sscal($s//R(1),$B->{$j}); } \%C }
# E = exp(Sum_k hbar^(2k) L_k):  E_m = (1/m) Sum_k k L_k E_{m-k}
my @E = ({ 0 => do { my $o=zero(); $o->[0]=R(1); $o } });
for my $m (1..$G){ my %acc; for my $k (1..$m){ my $t = pmul($L[$k],$E[$m-$k]); %acc = %{ padd(\%acc,$t,R($k)/$m) }; } $E[$m]=\%acc; }
# operators in t:  D = d/dX = t(1+t)/(1-tau t) d/dt;  dv/dX = -(1+t)/(t (1-tau t));  dv = -dt/t^2
my $one=zero(); $one->[0]=R(1);
my $pre = do { my $n=zero(); $n->[0]=R(1); $n->[1]=R(1); my $d=zero(); $d->[0]=R(1); $d->[1]=-$tau; # (1+t)/(1-tau t)
  my $di=zero(); for my $i (0..$T){ $di->[$i]=$tau**$i; } smul($n,$di) };
sub ddt { my ($f)=@_; my $r=zero(); for my $i (1..$T){ $r->[$i-1] = $i*$f->[$i]; } $r }
sub Dop { my ($f)=@_; my $g=ddt($f); my $r=zero(); for my $i (0..$T-1){ $r->[$i+1]=$g->[$i]; } smul($r,$pre) }   # t (1+t)/(1-tau t) f'
for my $g (1..$G){
  my $Eg = $E[$g];
  # F = int [u^1]E dv  + Sum_{r>=1} D^(r-1)([u^(r+1)]E dv/dX)
  my $F = zero();
  if (exists $Eg->{1}) { my $e=$Eg->{1}; # int -e(t) dt/t^2 : e = O(t^2g)
    for my $i (0..$T){ next if $e->[$i]->is_zero; my $pw=$i-2; die "log term" if $pw==-1; $F->[$pw+1] += -$e->[$i]/($pw+1) if $pw+1 <= $T && $pw+1>=0; } }
  for my $r (1..3*$g){ my $e=$Eg->{$r+1} or next; # e * dv/dX = -e (1+t)/(t(1-tau t))
    my $h = smul($e,$pre); my $f=zero(); for my $i (1..$T){ $f->[$i-1] = -$h->[$i]; } die "pole" unless $h->[0]->is_zero;
    for (1..$r-1){ $f = Dop($f); }
    $F = sadd($F,$f); }
  # [x^mu] F = (1/mu) Sum_k k f_k binom((tau+1)mu, mu-k)
  my @out; for my $mu (1..$M){ last if $mu > $T-2; my $c=R(0); for my $k (1..$mu){ next if $F->[$k]->is_zero; $c += $k*$F->[$k]*binomr(($tau+1)*$mu,$mu-$k); } $c /= $mu;
    my $xi = $c*($sig**($g-1))/binomr($mu*($tau+1)-1,$mu-1); push @out, "Xi_$g($mu)=".$xi; }
  print "g=$g: ", join("  ",@out), "\n";
  print "   F_t coefficients: ", join(" ", map { "t^$_:".$F->[$_] } grep { !$F->[$_]->is_zero } 0..($T>12?12:$T)), "\n" if $g==1;
}
