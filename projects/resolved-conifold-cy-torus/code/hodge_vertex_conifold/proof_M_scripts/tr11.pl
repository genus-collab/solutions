use strict; use warnings; use lib '.'; use ser; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
$ser::N = 10; my $N=$ser::N;
my $q = ser::zero(); for my $m (0..$N){ $q->[$m] = R(2*(-1)**$m, $m+2); }
my $w = ser::shiftup(ser::sqrt1($q),1);
my $T = ser::revert($w);
my $s = ser::compose($T, ser::scal(R(-1),$w));
print "s = ", ser::show($s,6), "\n";
my $sp = ser::deriv($s);
my $tt = ser::zero(); $tt->[1]=R(1);
my $A = ser::shiftdown(ser::sub_($tt,$s),1);
my $E = ser::zero(); for my $k (0..$N){ $E->[$k] = R(-(-1)**$k); }   # -1/(1+t)
for my $curve ('hurwitz','C') {
  my $Y = ser::zero();
  if ($curve eq 'hurwitz') { $Y->[1]=R(1); } else { for my $k (0..$N){ $Y->[$k] = R(-(-1)**$k); } }
  my $yd = ser::sub_($Y, ser::compose($Y,$s));
  my $Cd = ser::shiftdown($yd,1);
  my $den = ser::scal(R(2), ser::mul($Cd, ser::mul($E, ser::mul($A,$A))));
  my $H = ser::mul($sp, ser::inv($den));
  my @c;
  for my $k (1..5){ my $tk = ser::zero(); $tk->[$k]=R(1); my $sk = ser::zero(); $sk->[0]=R(1); for (1..$k){ $sk = ser::mul($sk,$s); }
    my $pr = ser::mul(ser::sub_($tk,$sk),$H); push @c, $pr->[3]; }
  print "$curve: omega11 = sum_k c_k du/u^(k+1), c_1..c_5 = ", join(", ", map { $_->bstr } @c), "\n";
  print "   times 24: ", join(", ", map { ($_*24)->bstr } @c), "\n";
}
