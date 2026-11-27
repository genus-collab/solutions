#!/usr/bin/perl
# Spot checks for the revised S3 text: Res_{z0} F_k (direct computation after (5)), formula (5), and the corner (17)
# against (16) with Theta_1 = -x^3/c^2. Exact rationals.
use strict; use warnings;
use Math::BigRat;
sub Q { return Math::BigRat->new("$_[0]"); }
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 1..$n; return $r; }
my ($pass, $fail) = (0, 0);
sub chk { my ($name, $a, $b) = @_; if ($a == $b) { $pass++; } else { $fail++; print "FAIL $name: $a vs $b\n"; } }
for my $xs ('3/7', '-2/5', '5/2') {
  my $x = Q($xs); my $c = $x + 1;
  for my $k (1..6) {
    # P_k = prod_{i=1}^k (k+1 - i c)^2
    my $P = Q(1); for my $i (1..$k) { my $f = Q($k+1) - $i*$c; $P *= $f*$f; }
    # direct: z0 = -c/(k+1); Res F_k = [prod (1+i z0)^2] / (k! z0^k prod_{i=2}^{k} (c + i z0)) * 1/(k+1)
    my $z0 = -$c/($k+1);
    my $num = Q(1); for my $i (1..$k) { my $f = 1 + $i*$z0; $num *= $f*$f; }
    my $den = fact($k) * ($z0 ** $k); for my $i (2..$k) { $den *= ($c + $i*$z0); }
    my $res = $num/$den/($k+1);
    my $resF = (Q(-1)**$k) * $P / ( Q($k+1)**2 * fact($k) * fact($k-1) * $c**(2*$k-1) );
    chk("Res F_k k=$k x=$xs", $res, $resF);
    # a_{k+1}^(k) = (k+1) (c z0/(c+z0)) Res F_k vs (5)
    my $a = ($k+1) * ($c*$z0/($c+$z0)) * $res;
    my $a5 = (Q(-1)**($k+1)) * $P / ( fact($k) * fact($k+1) * $c**(2*$k-2) );
    chk("(5) k=$k x=$xs", $a, $a5);
    my $g5 = (Q(-1)**$k) * $P / ( fact($k)**2 * $c**(2*$k-1) );
    chk("g top k=$k x=$xs", -($k+1)/$c * $a5, $g5);
  }
  # corner (17): e^(k)_{k1+1,k2+1} = (1/(k+2)) Theta_1 (k1+1)(k2+1) g_{k1+1}^(k1) g_{k2+1}^(k2)
  my $Th1 = -($x**3)/($c**2);
  for my $k (0..6) { for my $k1 (0..$k) { my $k2 = $k - $k1;
    my @g; for my $m ($k1, $k2) { my $P = Q(1); for my $i (1..$m) { my $f = Q($m+1) - $i*$c; $P *= $f*$f; }
      push @g, (Q(-1)**$m) * $P / ( fact($m)**2 * $c**(2*$m-1) ); }
    my $lhs = $Th1 * ($k1+1)*($k2+1) * $g[0]*$g[1] / ($k+2);
    my $P1 = Q(1); for my $i (1..$k1) { my $f = Q($k1+1) - $i*$c; $P1 *= $f*$f; }
    my $P2 = Q(1); for my $i (1..$k2) { my $f = Q($k2+1) - $i*$c; $P2 *= $f*$f; }
    my $rhs = (Q(-1)**($k+1)) * $x**3 * ($k1+1)*($k2+1) * $P1*$P2 / ( ($k+2) * fact($k1)**2 * fact($k2)**2 * $c**(2*$k) );
    chk("(17) k=$k k1=$k1 x=$xs", $lhs, $rhs);
  } }
  # leg t^1 lam_1 coefficient: (1-x)^2 x^3/c^(b+3) * (-1/2 + (2/3) 2^b - 1/4) = (1-x)^2 x^3 (2^(b+3)-9)/(12 c^(b+3))
  for my $b (0..7) {
    my $e11 = -((1-$x)**2)*$x**3/(2*$c**2); my $e12 = 2*((1-$x)**2)*$x**3/(3*$c**2);
    my $r1 = ((1-$x)**2)/(2*$c**2);
    my $lhs = ($e11 + (Q(2)**$b)*$e12)/($c**($b+1)) + $r1 * (-($x**3)/2)/($c**($b+1));
    my $rhs = ((1-$x)**2)*$x**3*((Q(2)**($b+3)) - 9)/(12*$c**($b+3));
    chk("leg t^1 b=$b x=$xs", $lhs, $rhs);
  }
}
print "passed $pass, failed $fail\n";
