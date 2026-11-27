package PJ;
# Univariate polynomials in j (array refs of Math::BigRat, index = power) and power series in t whose coefficients are
# such polynomials (array refs over the power of t).  Exact arithmetic.
use strict; use warnings; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
sub pz { return [Q(0)] }
sub pc { return [Q(shift)] }
sub pj { return [Q(0), Q(1)] }
sub ptrim { my $p = shift; pop @$p while @$p > 1 && $p->[-1] == 0; return $p }
sub padd { my ($a, $b, $cb) = @_; $cb //= Q(1); my $n = $#$a > $#$b ? $#$a : $#$b; return ptrim([ map { ($a->[$_] // Q(0)) + ($b->[$_] // Q(0)) * $cb } 0..$n ]) }
sub psub { my ($a, $b) = @_; return padd($a, $b, Q(-1)) }
sub pscale { my ($a, $c) = @_; return ptrim([ map { $_ * $c } @$a ]) }
sub pmul { my ($a, $b) = @_; my @c = map { Q(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i] == 0; for my $k (0..$#$b) { next if $b->[$k] == 0; $c[$i+$k] += $a->[$i] * $b->[$k] } } return ptrim(\@c) }
sub peval { my ($p, $x) = @_; my $s = Q(0); for (my $i = $#$p; $i >= 0; $i--) { $s = $s * $x + $p->[$i] } return $s }
sub pder { my $p = shift; return pz() if @$p <= 1; return ptrim([ map { $p->[$_] * $_ } 1..$#$p ]) }
sub iszero { my $p = shift; for (@$p) { return 0 if $_ != 0 } return 1 }
# falling factorial (j + c)(j + c - 1)...(j + c - n + 1) as a polynomial in j
sub ffj { my ($c, $n) = @_; my $p = pc(1); for my $i (0..$n-1) { $p = pmul($p, [Q($c - $i), Q(1)]) } return $p }
# its derivative in j
sub dffj { my ($c, $n) = @_; return pder(ffj($c, $n)) }
# exact division by (j - a); dies if not divisible
sub pdiv_lin { my ($p, $a) = @_; my $n = $#$p; return (pz(), $p->[0]) if $n == 0; my @q = map { Q(0) } 0..$n-1; my $acc = Q(0);
  for (my $i = $n; $i >= 0; $i--) { $acc = $acc * Q($a) + $p->[$i]; $q[$i-1] = $acc if $i >= 1; } return (ptrim(\@q), $acc) }
sub pstr { my $p = shift; my @t; for my $i (0..$#$p) { push @t, "$i:$p->[$i]" if $p->[$i] != 0 } return @t ? join(' ; ', @t) : '0' }
# ---- series in t (array over powers of t, entries are j-polynomials), truncated at order $N ----
our $N = 20;
sub sz { return [ map { pz() } 0..$N ] }
sub sone { my $s = sz(); $s->[0] = pc(1); return $s }
sub sadd { my ($a, $b, $c) = @_; return [ map { padd($a->[$_], $b->[$_], $c // Q(1)) } 0..$N ] }
sub sscalp { my ($a, $p) = @_; return [ map { pmul($a->[$_], $p) } 0..$N ] }   # multiply by a j-polynomial
sub smul { my ($a, $b) = @_; my $c = sz(); for my $i (0..$N) { next if iszero($a->[$i]); for my $k (0..$N-$i) { next if iszero($b->[$k]); $c->[$i+$k] = padd($c->[$i+$k], pmul($a->[$i], $b->[$k])) } } return $c }
sub sexp { my $a = shift; my $r = sz(); $r->[0] = pc(1); for my $n (1..$N) { my $s = pz(); for my $m (1..$n) { next if iszero($a->[$m]); $s = padd($s, pscale(pmul($a->[$m], $r->[$n-$m]), Q($m))) } $r->[$n] = pscale($s, Q(1)/$n) } return $r }
sub stder { my $a = shift; my $r = sz(); for my $n (0..$N) { $r->[$n] = pscale($a->[$n], Q($n)) } return $r }   # t d/dt
sub sderiv { my $a = shift; my $r = sz(); for my $n (0..$N-1) { $r->[$n] = pscale($a->[$n+1], Q($n+1)) } return $r }   # d/dt
sub sshift { my ($a, $k) = @_; my $r = sz(); for my $n (0..$N-$k) { $r->[$n+$k] = $a->[$n] } return $r }   # multiply by t^k
# scalar series (numbers) helpers
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
my @B = (Q(1)); sub bern { my $n = shift; while ($#B < $n) { my $m = scalar @B; my $s = Q(0); for my $k (0..$m-1) { $s += binom($m + 1, $k) * $B[$k] } push @B, -$s / ($m + 1) } return $B[$n] }
sub binom { my ($n, $k) = @_; return Q(0) if $k < 0 || $k > $n; my $r = Q(1); for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i } return $r }
# f = log S(t), S = (t/2)/sinh(t/2): f = Sum_k c_k t^{2k}, c_k = -B_{2k}/(2k (2k)!)
sub logS { my $s = sz(); for my $k (1..int($N/2)) { $s->[2*$k] = pc(-bern(2*$k) / (2*$k * fact(2*$k))) } return $s }
1;
