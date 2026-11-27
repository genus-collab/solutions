package SerS;
# Power series in L = lambda^2 whose coefficients are polynomials in sigma (array refs of BigRat, index = power of sigma).
# Basic objects of the r = 1 formula (6.1): e1 + e2 = 1, e1 e2 = -sigma, p_n = e1^n + e2^n,
#   log S(x) = Sum_k c_k x^{2k},  c_k = -B_{2k}/(2k (2k)!),
#   log(S1 S2) = Sum_k c_k p_{2k} L^k,   beta = e2 log S(e1 l) + e1 log S(e2 l) = -sigma Sum_k c_k p_{2k-1} L^k.
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
our $K = 8;     # truncation order in L
my @BN = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"), R("43867/798"), R("-174611/330"));
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 2..$n; $r }
sub ck { my $k = shift; return -$BN[$k] / (2*$k * fact(2*$k)); }
# sigma-polynomials
sub padd { my ($a, $b) = @_; my $n = $#$a > $#$b ? $#$a : $#$b; return [ map { ($a->[$_] // R(0)) + ($b->[$_] // R(0)) } 0..$n ]; }
sub pmul { my ($a, $b) = @_; my @c = map { R(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i] == 0; for my $j (0..$#$b) { $c[$i+$j] += $a->[$i] * $b->[$j] } } return \@c; }
sub pscale { my ($a, $c) = @_; return [ map { $_ * $c } @$a ]; }
my @P; sub pn { my $n = shift; if (!@P) { @P = ([R(2)], [R(1)]); } while ($#P < $n) { my $m = scalar @P; push @P, padd($P[$m-1], pmul([R(0), R(1)], $P[$m-2])); } return $P[$n]; }
# series: array ref (index = power of L) of sigma-polys
sub szero { return [ map { [R(0)] } 0..$K ]; }
sub sone { my $s = szero(); $s->[0] = [R(1)]; return $s; }
sub sadd { my ($a, $b) = @_; return [ map { padd($a->[$_], $b->[$_]) } 0..$K ]; }
sub sscale { my ($a, $c) = @_; return [ map { ref($c) eq 'ARRAY' ? pmul($a->[$_], $c) : pscale($a->[$_], $c) } 0..$K ]; }
sub smul { my ($a, $b) = @_; my $c = szero(); for my $i (0..$K) { for my $j (0..$K-$i) { $c->[$i+$j] = padd($c->[$i+$j], pmul($a->[$i], $b->[$j])); } } return $c; }
sub sexp { my $a = shift; # a[0] must be 0
  my $r = szero(); $r->[0] = [R(1)];
  for my $n (1..$K) { my $s = [R(0)]; for my $m (1..$n) { $s = padd($s, pscale(pmul($a->[$m], $r->[$n-$m]), R($m))); } $r->[$n] = pscale($s, R(1)/$n); }
  return $r; }
sub spow { my ($a, $n) = @_; my $r = sone(); $r = smul($r, $a) for 1..$n; return $r; }
sub logS1S2 { my $s = szero(); for my $k (1..$K) { $s->[$k] = pscale(pn(2*$k), ck($k)); } return $s; }
sub beta { my $s = szero(); for my $k (1..$K) { $s->[$k] = pmul([R(0), R(-1)], pscale(pn(2*$k-1), ck($k))); } return $s; }
sub S1S2 { return sexp(logS1S2()); }
# coefficient [L^h] as sigma-poly
sub coeff { my ($a, $h) = @_; return $a->[$h]; }
sub peval { my ($p, $s) = @_; my $r = R(0); my $x = R(1); for my $c (@$p) { $r += $c * $x; $x *= $s; } return $r; }
1;
