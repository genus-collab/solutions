use strict; use warnings; use Math::BigRat;
# Exact checks of the BS (v2) derivation behind eq:Rcf, Thm cfasym.
# (1) r_l from the (s1,s2) Feynman sum  [e^{z d1 d2/3} e^{s1^3/z} G](0,0),  G=(1-D'/27)^(-1/2),
#     D'(a,b)=a^2 b^2-4a^3+18ab-4b^3 (Disc at s3=1 is D'-27). BS claim: r0=1,r1=1/9,r2=2/81.
# (2) Stirling identities Theta(u,-1)=Theta(u)/(1-u), Theta(-u)Theta(u)=1, and
#     F2 assembly: Theta(-z/E,-1) / [Theta(-z/e1)Theta(-z/e2)] / h_w(z) = E/(E+z),
#     h_w = prod_{w in (E,-e1,-e2)} exp(-sum_n B_2n/(2n(2n-1)) w^{-(2n-1)} z^{2n-1})  (as in frozen.pl).
my $L = shift // 5; my $N = shift // 12;
sub R { Math::BigRat->new(@_) }
# ---- (1) bivariate series, weighted degree a:1, b:2, truncated at 3L
my $W = 3*$L;
sub bmul { my ($p,$q)=@_; my %r; for my $k1 (keys %$p){ my ($i1,$j1)=split /,/,$k1;
  for my $k2 (keys %$q){ my ($i2,$j2)=split /,/,$k2; my ($i,$j)=($i1+$i2,$j1+$j2); next if $i+2*$j>$W;
   $r{"$i,$j"} = ($r{"$i,$j"}//R(0)) + $p->{$k1}*$q->{$k2}; } } return \%r; }
my %Dp = ("2,2"=>R(1),"3,0"=>R(-4),"1,1"=>R(18),"0,3"=>R(-4));
my $x = { map { $_ => $Dp{$_}/27 } keys %Dp };   # x = D'/27, G=(1-x)^(-1/2)=sum binom(2n,n)/4^n x^n
my $G = {"0,0"=>R(1)}; my $xp = {"0,0"=>R(1)}; my $c = R(1);
for my $n (1..$W){ $xp = bmul($xp,$x); $c = $c*R(2*$n-1)/R(2*$n); last unless %$xp;
  for my $k (keys %$xp){ $G->{$k} = ($G->{$k}//R(0)) + $c*$xp->{$k}; } }
my @r = (R(0)) x ($L+1); my @fact=(R(1)); $fact[$_]=$fact[$_-1]*$_ for 1..3*$W;
for my $k (0..$W){ for my $m (0..int($k/3)){ my $l=$k-$m; next if $l>$L;
   my $co = $G->{($k-3*$m).",$k"} // R(0); next if $co==0;
   $r[$l] += R(1)/(R(3)**$k) * $fact[$k] * $co / $fact[$m]; } }
print "(1) r_l: ", join(", ", map {"r$_=".$r[$_]} 0..$L), "\n";
# ---- (2) Stirling series in z, exact, at given (e1,e2)
my @B=(R(1),R(-1,2),R(1,6),R(0),R(-1,30),R(0),R(1,42),R(0),R(-1,30),R(0),R(5,66),R(0),R(-691,2730),R(0),R(7,6),R(0),R(-3617,510));
sub binom { my ($n,$k)=@_; my $v=R(1); $v=$v*($n-$_+1)/$_ for 1..$k; $v }
sub bernpoly { my ($n,$x)=@_; my $s=R(0); for my $k (0..$n){ $s += binom($n,$k)*$B[$k]*($x**($n-$k)); } $s }
sub sexp { my @a=@_; my @e=(R(1),(R(0))x$N); for my $n(1..$N){ my $v=R(0); $v+=R($_)*$a[$_]*$e[$n-$_] for 1..$n; $e[$n]=$v/$n;} @e }
sub smul { my ($a,$b)=@_; my @c=(R(0))x($N+1); for my $i(0..$N){for my $j(0..$N-$i){$c[$i+$j]+=$a->[$i]*$b->[$j]}} \@c }
sub sinv { my ($a)=@_; my @b=(R(1)/$a->[0]); for my $n(1..$N){ my $v=R(0); $v+=$a->[$_]*$b[$n-$_] for 1..$n; $b[$n]=-$v/$a->[0]; } \@b }
# Theta(c z, x) = exp(-sum_k B_{k+1}(x) (-c z)^k /(k(k+1)))
sub theta { my ($cc,$x)=@_; my @l=(R(0)); for my $k(1..$N){ $l[$k] = -bernpoly($k+1,$x)*((-$cc)**$k)/R($k*($k+1)); } [sexp(@l)] }
for my $pr ([2,1],[5,11],[7,3],[1,-3]) { my ($e1,$e2)=map{R($_)}@$pr; my $E=$e1+$e2;
  # identity checks at u = c z with c=-1/E
  my $c=R(-1)/$E; my $t1=theta($c,R(-1)); my $t0=theta($c,R(0));
  my $oneminus=[R(1),-$c,(R(0))x($N-1)]; my $chk=smul($t1,$oneminus); my $ok1=1; $ok1&&=($chk->[$_]==$t0->[$_]) for 0..$N;
  my $ok2=1; my $pp=smul(theta($c,R(0)),theta(-$c,R(0))); $ok2&&=($pp->[$_]==($_==0?1:0)) for 0..$N;
  # h_w as in frozen.pl
  my @lh=(R(0))x($N+1); for my $n(1..int(($N+1)/2)){ my $m=2*$n-1; next if $m>$N; my $S=R(0); $S+= R(1)/($_**$m) for ($E,-$e1,-$e2); $lh[$m]=-$B[2*$n]/R(2*$n*(2*$n-1))*$S; }
  my $h=[sexp(@lh)];
  my $Ups = smul(smul(theta(R(1)/$e1,R(0)),theta(R(1)/$e2,R(0))), sinv(theta(R(1)/$E,R(0))));  # Upsilon(z/e1,z/e2)
  my $ok3=1; $ok3&&=($Ups->[$_]==$h->[$_]) for 0..$N;
  my $num=theta($c,R(-1)); my $den=smul(smul(theta(R(-1)/$e1,R(0)),theta(R(-1)/$e2,R(0))),$h);
  my $res=smul($num,sinv($den)); my $ok4=1; $ok4&&=($res->[$_]==(R(-1)/$E)**$_) for 0..$N;
  print "(2) eps=(@$pr): Theta(u,-1)(1-u)=Theta(u): ".($ok1?"OK":"FAIL").";  Theta(u)Theta(-u)=1: ".($ok2?"OK":"FAIL").
        ";  h_w==Upsilon(z/e1,z/e2): ".($ok3?"OK":"FAIL")."; Theta(-z/E,-1)/(N0 h_w) == E/(E+z) through z^$N: ".($ok4?"OK":"FAIL")."\n";
}
