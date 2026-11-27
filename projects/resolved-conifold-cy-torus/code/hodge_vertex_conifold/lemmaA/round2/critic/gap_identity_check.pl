#!/usr/bin/perl
# Critic check: R2's Corollary 7 (sigma^h Int_{Mbar_h} K Omega = -(2h-3)! [l^{2h}] S(e1 l) S(e2 l), c = 1, e1 = -tau, e2 = 1+tau)
# converted to the frame-B normalisation c = 1 + x, tau = -x/c (Xi^c_h(empty) = c^-(3h-3) Xi_h(empty)), gives
#     (-x c^3)^(h-1) Z[h,0] = (2h-3)! P_h(x,1)/x = gap_h,   P_h(e1,e2) = [l^{2h}] S(e1 l) S(e2 l)   (rc_lib.m Pg, cyPolar),
# i.e. F4 at n = 0 in every genus.  This script checks the algebra of that conversion exactly (no data involved).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $HMAX = shift // 10;
my $N = 2*$HMAX;
# S(y) = (y/2)/sinh(y/2) as a power series in y (even); 1/S = sum_k (y/2)^(2k)/(2k+1)!
my @inv = map { R(0) } 0..$N;
my $f = R(1);
for my $k (0..$HMAX) { my $fact = R(1); $fact *= $_ for 1..(2*$k+1); $inv[2*$k] = R(1)/($fact * R(4)**$k); }
my @S = map { R(0) } 0..$N; $S[0] = R(1);
for my $n (1..$N) { my $s = R(0); for my $k (1..$n) { $s += $inv[$k]*$S[$n-$k]; } $S[$n] = -$s; }
sub Pcoef {   # [l^{2h}] S(a l) S(b l)
  my ($h,$a,$b) = @_; my $s = R(0);
  for my $i (0..2*$h) { $s += $S[$i]*$a**$i * $S[2*$h-$i]*$b**(2*$h-$i); } $s }
sub fact { my $n=shift; my $f=R(1); $f*=$_ for 1..$n; $f }
my $bad = 0; my $cnt = 0;
for my $x (R(2), R("1/3"), R("-3/5"), R("7/4")) {
  my $c = 1 + $x; my $tau = -$x/$c; my $sig = $tau*(1+$tau); my ($e1,$e2) = (-$tau, 1+$tau);
  for my $h (2..$HMAX) {
    my $IKO = -fact(2*$h-3) * Pcoef($h,$e1,$e2) / $sig**$h;          # Int K Omega at c = 1 (R2 Cor 7)
    my $Z   = $IKO / $c**(3*$h-3);                                    # frame-B normalisation
    my $lhs = (-$x*$c**3)**($h-1) * $Z;
    my $gap = fact(2*$h-3) * Pcoef($h,$x,R(1)) / $x;                  # [t^{2-2h}] P_h Li_{3-2h}(e^t)/(e1 e2), e1 = x, e2 = 1
    $cnt++; if ($lhs != $gap) { $bad++; print "MISMATCH x=$x h=$h lhs=$lhs gap=$gap\n"; }
  }
}
print "gap identity: $cnt comparisons, $bad mismatches\n";
