# norlund_check.pl: checks the two Norlund-type identities used in the proof of (6.1) from the closed one-point formula,
# at NON-integer N = j tau (where they are proved by polynomiality):
#   A(t) := exp( Sum_k t^(2k) A_k ),  A_k = (2k-1)! (j/al0) c_k(j,al0),  al0 = 1/tau      =  Sum_n (-N)_n [y^n] S1(y)^(1-N) (t/tau)^n
#   B(t) := exp( Sum_k t^(2k) B_k ),  B_k = (2k-1)! (j/al1) c_k(j,al1),  al1 = -1/(tau+1) =  Sum_n ( M)_n [y^n] S1(y)^(1+M) (t/(tau+1))^n
# with c_k(u,al) = [t^(2k)] S(ut)/S(al t), S(t) = 2 sinh(t/2)/t, S1 = 1/S, N = j tau, M = j(tau+1), (a)_n = a(a-1)...(a-n+1).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $tau = R(shift // "3/7"); my $JM = shift // 6; my $KM = shift // 6;
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
my @s = map { R(1)/(R(4)**$_ * fact(2*$_+1)) } 0..$KM+1;
my @sp=(R(1)); for my $n (1..$KM+1){ my $acc=R(0); $acc += $s[$_]*$sp[$n-$_] for 1..$n; push @sp, -$acc; }
sub ck { my ($u,$a,$k)=@_; my $c=R(0); for my $p (0..$k){ $c += $s[$p]*$sp[$k-$p]*($u**(2*$p))*($a**(2*($k-$p))); } $c }
# [y^(2m)] S1(y)^e : log S1 = Sum_k b_k y^(2k)
my @B=(R(1)); for my $n (1..2*$KM+2){ my $acc=R(0); for my $k (0..$n-1){ my $bin=fact($n+1)/(fact($k)*fact($n+1-$k)); $acc += $bin*$B[$k]; } push @B, -$acc/($n+1); }
my @b=(R(0)); for my $k (1..$KM+1){ push @b, -$B[2*$k]/(2*$k*fact(2*$k)); }
sub S1pow { my ($e,$m)=@_; my @ex=(R(1)); for my $q (1..$m){ my $acc=R(0); $acc += $_*$e*$b[$_]*$ex[$q-$_] for 1..$q; push @ex, $acc/$q; } $ex[$m] }
sub falling { my ($a,$n)=@_; my $r=R(1); $r *= ($a-$_) for 0..$n-1; $r }
sub expser { my ($c)=@_; my @ex=(R(1)); for my $q (1..$KM){ my $acc=R(0); $acc += $_*$c->[$_]*$ex[$q-$_] for 1..$q; push @ex, $acc/$q; } \@ex }   # exp(Sum c_k t^2k)
my ($n,$bad)=(0,0);
for my $j (1..$JM){ my $N=$j*$tau; my $M=$j*($tau+1); my @al=(R(1)/$tau, R(-1)/($tau+1));
  my @Ak=(R(0)); my @Bk=(R(0)); for my $k (1..$KM){ push @Ak, fact(2*$k-1)*($j/$al[0])*ck(R($j),$al[0],$k); push @Bk, fact(2*$k-1)*($j/$al[1])*ck(R($j),$al[1],$k); }
  my $EA=expser(\@Ak); my $EB=expser(\@Bk);
  for my $m (0..$KM){ my $a = falling(-$N,2*$m)*S1pow(1-$N,$m)/($tau**(2*$m)); my $bb = falling($M,2*$m)*S1pow(1+$M,$m)/(($tau+1)**(2*$m)); $n+=2;
    if ($a != $EA->[$m]) { $bad++; print "A mismatch j=$j m=$m: ".$EA->[$m]." vs $a\n"; }
    if ($bb != $EB->[$m]) { $bad++; print "B mismatch j=$j m=$m: ".$EB->[$m]." vs $bb\n"; } } }
print "tau=$tau: $n Norlund-identity checks (j<=$JM, m<=$KM), $bad mismatches\n";
