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
  # alg = (delta/c)^(1/2) / (c s*)
  my $alg = mul(powS($dc,R(1,2)), inv(scl($c,$s)));
  # alg^2 = delta/(c^3 s*^2)
  chk("alg^2 x=$xs", mul($alg,$alg), mul(scl(R(1)/($c**3),$delta), inv(mul($s,$s))));
  # branch-free alg = (1-u)^(1/4) (1 + (1-u)^(1/2))/2
  my $oneu = sub_(cst(1), scl(4*$x/($c*$c), $tau));
  my $algbf = scl(R(1,2), mul(powS($oneu,R(1,4)), add(cst(1), powS($oneu,R(1,2)))));
  chk("alg branch-free x=$xs", $alg, $algbf);
  # rho = tau/(t alg); check rho_1 and rho_2
  my $tau_t = divt($tau);
  my $rho = mul($tau_t, inv($alg));
  my $r1 = -((1-$x)**2)/(2*$c*$c);
  my $r2 = (1 - 8*$x + 15*$x**2 - 8*$x**3 + $x**4)/(6*$c**4);
  chk("rho_1,rho_2 x=$xs", $rho, [R(1),$r1,$r2,(map {$rho->[$_]} 3..$N)], 2);
  my $rinv=inv($rho);
  chk("rho^-1_2 x=$xs", $rinv, [R(1), (1-$x)**2/(2*$c*$c), (1+4*$x-12*$x**2+4*$x**3+$x**4)/(12*$c**4), (map {$rinv->[$_]} 3..$N)], 2);
  # L_1, L_2, Delta_1, Delta_2 from the literal formulas: L_1 = c/tau - x s*, L_2 = x s*
  my $L2 = scl($x,$s);
  # t*L_1 = c t/tau - x t s*
  my $tL1 = sub_(scl($c, inv($tau_t)), scl($x, mul(tvar(),$s)));
  # Delta_2 = x^3 delta s*^2 ; also Delta_2 = (L_2 - L_1)(x^2 - c x L_2): times t: t(L_2-L_1) = t L_2 - t L_1
  my $D2 = scl($x**3, mul($delta, mul($s,$s)));
  my $tL2mL1 = sub_(mul(tvar(),$L2), $tL1);
  my $x2mcxL2 = sub_(cst($x*$x), scl($c*$x,$L2));
  # x^2 - c x L_2 = -x^3 tau s*^2  (O(t)); so (L_2-L_1)(x^2 - c x L_2) = t(L2-L1) * [(x^2-cxL2)/t]
  chk("x^2-cxL2 = -x^3 tau s*^2 x=$xs", $x2mcxL2, scl(-$x**3, mul($tau, mul($s,$s))));
  chk("Delta_2 = (L2-L1)(x^2-cxL2) x=$xs", $D2, mul($tL2mL1, divt($x2mcxL2)));
  # Delta_2 expansion: x^3/c - x^5 tau^2/c^5 + O(tau^3)
  my $D2tau = add(cst($x**3/$c), scl(-$x**5/$c**5, mul($tau,$tau)));
  chk("Delta_2 to tau^2 x=$xs", $D2, $D2tau, 2);
  # t^2 rho^2 Delta_1 = -x c^3 with Delta_1 = -x delta/(tau^2 s*^2): t^2 Delta_1 = -x delta/((tau/t)^2 s*^2)
  my $t2D1 = scl(-$x, mul($delta, inv(mul(mul($tau_t,$tau_t), mul($s,$s)))));
  chk("t^2 rho^2 Delta_1 = -x c^3 x=$xs", mul(mul($rho,$rho),$t2D1), cst(-$x*$c**3));
  # Delta_1 = (L_1 - L_2)(x^2 - c x L_1): t^2 Delta_1 = t(L1-L2) * t(x^2 - c x L1)
  my $tx2mcxL1 = sub_(scl($x*$x, tvar()), scl($c*$x, $tL1));
  chk("Delta_1 = (L1-L2)(x^2-cxL1) x=$xs", $t2D1, mul(scl(R(-1),$tL2mL1), $tx2mcxL1));
  # Theta = (L_2 - L_1)(c x L_2 - x^2)^2/Delta_2 = x^2 - c x L_2 = x c dphi' - x ; dphi' = 1 - x s*
  my $Theta_def = mul(divt(mul($tL2mL1, mul($x2mcxL2,$x2mcxL2))), inv($D2));
  # careful: (L2-L1) = tL2mL1/t ; times (x^2-cxL2)^2 which is O(t^2): divide by t once
  chk("Theta def = x^2 - c x L2 x=$xs", $Theta_def, $x2mcxL2);
  my $dphip = sub_(cst(1), scl($x,$s));
  chk("Theta = x c dphi' - x x=$xs", $x2mcxL2, sub_(scl($x*$c,$dphip), cst($x)));
  # dphi from definition: phi_tau(s*) - phi_0(1/c) = log(c s*) + x log((1-s*) c/x) - log(1 - tau s*)
  my $dphi_def = add(add(logS(scl($c,$s)), scl($x, logS(scl($c/$x, sub_(cst(1),$s))))), scl(R(-1), logS(sub_(cst(1), mul($tau,$s)))));
  chk("dphi = Int (1 - x s*) dt x=$xs", $dphi_def, integ($dphip));
  # e^dphi = (c s*)^c U_1^(x-1): dphi = c log(c s*) + (x-1) log(1 - tau s*)
  chk("dphi = c log(c s*) + (x-1) log U1 x=$xs", $dphi_def, add(scl($c, logS(scl($c,$s))), scl($x-1, logS(sub_(cst(1),mul($tau,$s))))));
  # dphi low orders in t: t/c - x^2 t^2/(2c^3) + (1-x)^2 x^2 t^3/(6 c^5)
  chk("dphi low orders x=$xs", $dphi_def, [R(0), R(1)/$c, -$x**2/(2*$c**3), (1-$x)**2*$x**2/(6*$c**5), (map {$dphi_def->[$_]} 4..$N)], 3);
  # l = c(1/tau - 1/t) - x s*  ; t l = c t/tau - c - x t s* ; l_0 = c/2 - x/c = (1+x^2)/(2c), l_1 = c/12 - x^2/c^3
  my $tl = sub_($tL1, cst($c));
  my $l = divt($tl);
  chk("l_0, l_1 x=$xs", $l, [(1+$x*$x)/(2*$c), $c/12 - $x**2/$c**3, (map {$l->[$_]} 2..$N)], 1);
  # zeta = tau e^{dphi/c} = c (1-U1) U1^{-(1-x)/c}
  my $zeta = mul($tau, expS(scl(R(1)/$c,$dphi_def)));
  my $U1 = sub_(cst(1), mul($tau,$s));
  my $zeta2 = scl($c, mul(mul($tau,$s), powS($U1, -(1-$x)/$c)));
  chk("zeta two forms x=$xs", $zeta, $zeta2);
  chk("zeta low x=$xs", $zeta, [R(0), R(1), R(1)/($c*$c) - R(1,2), (map {$zeta->[$_]} 3..$N)], 2);
  # Theta series: -x^3 t/c^2 + x^3 (1-x)^2 t^2/(2 c^4) - x^3 (1 - 8x + 12x^2 - 8x^3 + x^4) t^3/(6 c^6)
  chk("Theta series x=$xs", $x2mcxL2, [R(0), -$x**3/$c**2, $x**3*(1-$x)**2/(2*$c**4), -$x**3*(1-8*$x+12*$x**2-8*$x**3+$x**4)/(6*$c**6), (map {$x2mcxL2->[$_]} 4..$N)], 3);
}
print "spot checks: $nchk, failures: $fail\n";
