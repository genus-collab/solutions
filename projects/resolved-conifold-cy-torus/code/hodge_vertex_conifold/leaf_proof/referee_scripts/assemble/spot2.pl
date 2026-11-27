#!/usr/bin/perl
# Exact spot checks (Math::BigRat) of formulas restated in leaf_structure_proof.md in unified notation.
# Truncated power series in t with rational coefficients at rational x.
use strict; use warnings;
use Math::BigRat;
my $N = 9;   # truncation order in t (keep t^0..t^N)
sub R { Math::BigRat->new(@_) }
sub zero { [map { R(0) } 0..$N] }
sub add { my ($a,$b)=@_; [map { $a->[$_] + $b->[$_] } 0..$N] }
sub sub_ { my ($a,$b)=@_; [map { $a->[$_] - $b->[$_] } 0..$N] }
sub scl { my ($k,$a)=@_; [map { $k * $a->[$_] } 0..$N] }
sub mul { my ($a,$b)=@_; my $c=zero(); for my $i (0..$N){ next if $a->[$i]==0; for my $j (0..$N-$i){ $c->[$i+$j] += $a->[$i]*$b->[$j]; } } $c }
sub inv { my ($a)=@_; die "inv: zero constant" if $a->[0]==0; my $c=zero(); $c->[0]=R(1)/$a->[0];
  for my $n (1..$N){ my $s=R(0); for my $k (1..$n){ $s += $a->[$k]*$c->[$n-$k]; } $c->[$n] = -$s/$a->[0]; } $c }
sub cst { my ($v)=@_; my $c=zero(); $c->[0]=R($v); $c }
sub tvar { my $c=zero(); $c->[1]=R(1); $c }
sub der { my ($a)=@_; my $c=zero(); for my $n (0..$N-1){ $c->[$n]=($n+1)*$a->[$n+1]; } $c }
sub integ { my ($a)=@_; my $c=zero(); for my $n (0..$N-1){ $c->[$n+1]=$a->[$n]/($n+1); } $c }
sub expS { my ($a)=@_; die "exp: nonzero constant" unless $a->[0]==0; my $c=cst(1); my $term=cst(1);
  for my $k (1..$N){ $term = scl(R(1)/$k, mul($term,$a)); $c=add($c,$term); } $c }
sub logS { my ($a)=@_; die "log: constant not 1" unless $a->[0]==1; integ(mul(der($a),inv($a))) }
sub powS { my ($a,$e)=@_; # a with constant 1, rational exponent e
  my $c=cst(1); my $u=sub_($a,cst(1)); my $term=cst(1); my $coef=R(1);
  for my $k (1..$N){ $coef = $coef*($e-($k-1))/$k; $term=mul($term,$u); $c=add($c,scl($coef,$term)); } $c }
sub sqrtS { my ($a)=@_; die "sqrt: constant not 1" unless $a->[0]==1; powS($a, R(1,2)) }
sub divt { my ($a)=@_; die "divt: constant nonzero" unless $a->[0]==0; my $c=zero(); for my $n (0..$N-1){ $c->[$n]=$a->[$n+1]; } $c->[$N]=R(0); $c }
my $fail=0; my $nchk=0;
sub chk { my ($name,$a,$b,$upto)=@_; $upto //= $N-2; $nchk++; for my $n (0..$upto){ if ($a->[$n] != $b->[$n]) { print "FAIL $name at t^$n: ",$a->[$n]," vs ",$b->[$n],"\n"; $fail++; return; } } }

for my $xs ('2','3/5','-5/2','1/3','7/3') {
  my $x=R($xs); my $c=$x+1;
  # tau = 1 - e^{-t}
  my $tau = sub_(cst(1), expS(scl(R(-1),tvar())));
  # s* from fixed point s = (1 + x tau s^2)/c
  my $s=cst(R(1)/$c); for (1..$N+2){ $s = scl(R(1)/$c, add(cst(1), scl($x, mul($tau, mul($s,$s))))); }
  chk("critical eq x=$xs", add(sub_(scl($x,mul($tau,mul($s,$s))), scl($c,$s)), cst(1)), zero());
  # delta/c = sqrt(1 - 4 x tau/c^2)
  my $dc = sqrtS(sub_(cst(1), scl(4*$x/($c*$c), $tau)));
  my $delta = scl($c,$dc);
  chk("s* = 2/(c+delta) x=$xs", $s, scl(R(2), inv(add(cst($c),$delta))));
for my $xs ('2','3/5','-5/2','1/3','7/3') {
  my $x=R($xs); my $c=$x+1;
  my $tau = sub_(cst(1), expS(scl(R(-1),tvar())));
  my $s=cst(R(1)/$c); for (1..$N+2){ $s = scl(R(1)/$c, add(cst(1), scl($x, mul($tau, mul($s,$s))))); }
  my $dc = sqrtS(sub_(cst(1), scl(4*$x/($c*$c), $tau)));
  my $alg = mul(powS($dc,R(1,2)), inv(scl($c,$s)));
  my $rho = mul(divt($tau), inv($alg));
  my $r3 = -((1-$x)**2)*(1 - 20*$x + 60*$x**2 - 20*$x**3 + $x**4)/(24*$c**6);
  my $r4 = (1 - 52*$x + 493*$x**2 - 1624*$x**3 + 2365*$x**4 - 1624*$x**5 + 493*$x**6 - 52*$x**7 + $x**8)/(120*$c**8);
  chk("rho_3,rho_4 x=$xs", $rho, [$rho->[0],$rho->[1],$rho->[2],$r3,$r4,(map {$rho->[$_]} 5..$N)], 4);
  my $rm2 = inv(mul($rho,$rho));
  my $d2 = (5 - 4*$x - 6*$x**2 - 4*$x**3 + 5*$x**4)/(12*$c**4);
  my $d3 = (1-$x)**2/(12*$c**2);
  chk("Delta_1 bracket x=$xs", $rm2, [R(1), (1-$x)**2/$c**2, $d2, $d3, (map {$rm2->[$_]} 4..$N)], 3);
  # dphi t^4 coefficient: -x^2 (1 - 8x + 12x^2 - 8x^3 + x^4)/(24 c^7)
  my $dphip = sub_(cst(1), scl($x,$s));
  my $dphi = integ($dphip);
  chk("dphi_4 x=$xs", $dphi, [(map {$dphi->[$_]} 0..3), -$x**2*(1-8*$x+12*$x**2-8*$x**3+$x**4)/(24*$c**7), (map {$dphi->[$_]} 5..$N)], 4);
  # A_j(t) = (-1)^j C_j t^-1 zeta^j 2F1(jx/c, jx/c; j; tau); check j=1 at t^1, t^2: a_1^(1) = -(1-x)^2/2, a_1^(2) = -(2-4x+3x^2-4x^3+2x^4)/(12c^2)
  my $zeta = mul($tau, expS(scl(R(1)/$c, $dphi)));
  for my $j (1..4) {
    my $Cj = R(1); for my $i (1..$j-1){ $Cj *= ($j - $i*$c)**2; } my $fact=R(1); for my $i (1..$j){ $fact*=$i; } my $fact1=R(1); for my $i (1..$j-1){ $fact1*=$i; }
    $Cj = $Cj/($fact*$fact1*($c**(2*$j-4)));
    # 2F1(a,a;j;tau) series
    my $a = $j*$x/$c; my $F=cst(0); my $term=cst(1); my $coef=R(1);
    my $taun=cst(1);
    for my $n (0..$N){ $F = add($F, scl($coef,$taun)); $coef = $coef*($a+$n)*($a+$n)/(($j+$n)*($n+1)); $taun=mul($taun,$tau); }
    my $zj=cst(1); for (1..$j){ $zj=mul($zj,$zeta); }
    my $Aj = scl(((-1)**$j)*$Cj, divt(mul($zj,$F)));
    if ($j==1) { chk("A_1 low x=$xs", $Aj, [-$c*$c, -((1-$x)**2)/2, -(2-4*$x+3*$x**2-4*$x**3+2*$x**4)/(12*$c*$c), (map {$Aj->[$_]} 3..$N)], 2); }
    if ($j==2) { chk("A_2 low x=$xs", $Aj, [R(0), (1-$x)**2/2, (1-$x)**4/(2*$c*$c), (map {$Aj->[$_]} 3..$N)], 2); }
    if ($j==3) { chk("A_3 low x=$xs", $Aj, [R(0), R(0), -(2-$x)**2*(1-2*$x)**2/(12*$c*$c), -(2-$x)**2*(1-2*$x)**2*(1-$x)**2/(8*$c**4), (map {$Aj->[$_]} 4..$N)], 3); }
    if ($j==4) { chk("A_4 top x=$xs", $Aj, [R(0), R(0), R(0), (3-$x)**2*(1-$x)**2*(1-3*$x)**2/(36*$c**4), (map {$Aj->[$_]} 4..$N)], 3); }
  }
}
print "spot2 checks: $nchk, failures: $fail\n";
