use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
my @B = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"));  # B_0,B_2,..,B_10 (index k -> B_{2k})
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
my $K = 5; # up to lambda^{2K}
# log S(x) = - sum_k B_{2k} x^{2k} / (2k (2k)!)
my @lS = (R(0)); for my $k (1..$K){ push @lS, -$B[$k]/(2*$k*fact(2*$k)); }
sub pexp { # exp of series in u=lambda^2 with zero constant term, coefficients c[1..K]
  my @c=@_; my @r=(R(1)); push @r, R(0) for 1..$K;
  # r' = c' r  => n r_n = sum_{m=1}^n m c_m r_{n-m}
  for my $n (1..$K){ my $s=R(0); for my $m (1..$n){ $s += $m*$c[$m]*$r[$n-$m]; } $r[$n]=$s/$n; }
  return @r;
}
sub formula { my ($e1,$e2,$j,$h)=@_;
  my @c=(R(0)); for my $k (1..$K){ push @c, (1+$j*$e2)*$lS[$k]*($e1**(2*$k)) + (1+$j*$e1)*$lS[$k]*($e2**(2*$k)); }
  my @r=pexp(@c); return $r[$h]; }
sub data { my ($sig,$j,$h)=@_;
  if($h==1){ return (-($sig*($j-2)-1))/24 * (-1) ; } # placeholder replaced below
}
sub psihat { my ($sig,$j,$h)=@_;
  if($h==1){ # Xi_1 = (j-2-1/sigma)/24, psihat_1 = sigma*Xi_1
     return $sig*($j-2-1/$sig)/24; }
  if($h==2){ return (7 + $sig*(28-12*$j) + $sig**2*(5*$j**2-26*$j+24))/5760; }
  if($h==3){ return (-93 + $sig*(205*$j-558) - $sig**2*(147*$j**2-878*$j+984) + $sig**3*(35*$j**3-336*$j**2+836*$j-480))/2903040; }
  if($h==4){ return (1143 + $sig*(9144-2912*$j) + $sig**2*(24720-18524*$j+2714*$j**2) + $sig**3*(25728-35568*$j+12184*$j**2-1120*$j**3) + $sig**4*(8064-19856*$j+12436*$j**2-2660*$j**3+175*$j**4))/1393459200; }
}
my $bad=0; my $n=0;
for my $tau (R("-1/3"), R("2/5"), R("3"), R("-7/2"), R("1/7")) {
  my $e1 = -$tau; my $e2 = 1+$tau; my $sig = $tau*(1+$tau);
  for my $j (0..6, R("1/2"), R("-3")) { $j = R($j) unless ref $j;
    for my $h (1..4) { $n++;
      my $f = formula($e1,$e2,$j,$h); my $d = psihat($sig,$j,$h);
      if ($f != $d) { $bad++; print "MISMATCH tau=$tau j=$j h=$h formula=$f data=$d\n"; } } } }
print "tested $n points, mismatches $bad\n";
