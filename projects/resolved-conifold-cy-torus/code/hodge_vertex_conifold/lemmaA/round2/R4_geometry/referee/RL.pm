package RL;
# Referee library (written from scratch for the R4 referee pass; shares no code with ../PJ.pm, ../RefHodge.pm).
# Exact arithmetic with Math::BigRat.
#  - univariate polynomials in j: array ref, index = power of j
#  - bivariate polynomials in (j, y): array ref indexed by the power of y, entries univariate polynomials in j
#  - scalar power series in t: array ref, index = power of t
use strict; use warnings; use Math::BigRat; use Math::BigInt;
sub Q { return Math::BigRat->new($_[0]) }
my $ZERO = Q(0);
# ---------- numbers ----------
my @FACT = (Q(1));
sub fact { my $n = shift; while (@FACT <= $n) { push @FACT, $FACT[-1] * Q(scalar @FACT) } return $FACT[$n] }
sub binq { my ($x, $k) = @_; return Q(0) if $k < 0; my $r = Q(1); for my $i (0..$k-1) { $r *= (Q($x) - $i) } return $r / fact($k) }
my @BN = (Q(1));
# Bernoulli numbers from the generating function t/(e^t - 1): Sum_{k=0}^{n} binom(n+1,k) B_k = 0 for n >= 1
sub bern { my $n = shift; while (@BN <= $n) { my $m = scalar @BN; my $s = Q(0); for my $k (0..$m-1) { $s += binq($m + 1, $k) * $BN[$k] } push @BN, -$s / Q($m + 1) } return $BN[$n] }
# ---------- univariate polynomials ----------
sub ptrim { my $p = shift; pop @$p while @$p > 1 && $p->[-1] == 0; $p = [Q(0)] unless @$p; return $p }
sub pc { return [Q($_[0])] }
sub padd { my ($a, $b, $c) = @_; my $n = @$a > @$b ? scalar @$a : scalar @$b; my @r;
  for my $i (0..$n-1) { my $x = defined $a->[$i] ? $a->[$i] : $ZERO; my $y = defined $b->[$i] ? $b->[$i] : $ZERO; push @r, defined $c ? $x + $y * $c : $x + $y } return ptrim(\@r) }
sub psc { my ($a, $c) = @_; return ptrim([ map { $_ * $c } @$a ]) }
sub pmul { my ($a, $b) = @_; my @r = map { Q(0) } 0..($#$a + $#$b);
  for my $i (0..$#$a) { next if $a->[$i] == 0; for my $k (0..$#$b) { next if $b->[$k] == 0; $r[$i+$k] += $a->[$i] * $b->[$k] } } return ptrim(\@r) }
sub pev { my ($p, $x) = @_; my $s = Q(0); for (my $i = $#$p; $i >= 0; $i--) { $s = $s * $x + $p->[$i] } return $s }
sub pder { my $p = shift; return [Q(0)] if @$p <= 1; return ptrim([ map { $p->[$_] * $_ } 1..$#$p ]) }
sub pzero { my $p = shift; for (@$p) { return 0 if $_ != 0 } return 1 }
sub peq { my ($a, $b) = @_; return pzero(padd($a, $b, Q(-1))) }
# (j + c)_(n) = prod_{i=0}^{n-1} (j + c - i); for n < 0 not supported
sub ffj { my ($c, $n) = @_; my $p = [Q(1)]; for my $i (0..$n-1) { $p = pmul($p, [Q($c) - $i, Q(1)]) } return $p }
# divide by (j - a): returns (quotient, remainder)
sub pdivlin { my ($p, $a) = @_; my $n = $#$p; return ([Q(0)], $p->[0]) if $n == 0; my @q = map { Q(0) } 0..$n-1; my $acc = Q(0);
  for (my $i = $n; $i >= 1; $i--) { $acc = $acc * $a + $p->[$i]; $q[$i-1] = $acc } my $rem = $acc * $a + $p->[0]; return (ptrim(\@q), $rem) }
sub pstr { my $p = shift; my @t; for my $i (0..$#$p) { push @t, "$i:" . $p->[$i] if $p->[$i] != 0 } return @t ? join(' ; ', @t) : '0' }
# ---------- bivariate polynomials P(j,y) = Sum_b P[b](j) y^b ----------
sub btrim { my $p = shift; pop @$p while @$p > 1 && pzero($p->[-1]); return $p }
sub badd { my ($a, $b, $c) = @_; my $n = @$a > @$b ? scalar @$a : scalar @$b; my @r;
  for my $i (0..$n-1) { push @r, padd($a->[$i] // [Q(0)], $b->[$i] // [Q(0)], $c) } return btrim(\@r) }
sub bsc { my ($a, $c) = @_; return btrim([ map { psc($_, $c) } @$a ]) }
sub bmul { my ($a, $b) = @_; my @r = map { [Q(0)] } 0..($#$a + $#$b);
  for my $i (0..$#$a) { next if pzero($a->[$i]); for my $k (0..$#$b) { next if pzero($b->[$k]); $r[$i+$k] = padd($r[$i+$k], pmul($a->[$i], $b->[$k])) } } return btrim(\@r) }
sub bmulj { my ($a, $p) = @_; return btrim([ map { pmul($_, $p) } @$a ]) }   # multiply by a polynomial in j
# (j + y + c)_(n)
sub ffjy { my ($c, $n) = @_; my $p = [[Q(1)]]; for my $i (0..$n-1) { $p = bmul($p, [[Q($c) - $i, Q(1)], [Q(1)]]) } return $p }
# (y + c)_(n) as bivariate (constant in j)
sub ffy { my ($c, $n) = @_; my $u = ffj($c, $n); return btrim([ map { [$_] } @$u ]) }
# divide P(j,y) by (j + y - a), i.e. by (y - r(j)) with r(j) = a - j; returns (quotient, remainder(j))
sub bdivlin { my ($P, $a) = @_; my $n = $#$P; my $r = [Q($a), Q(-1)]; return ([[Q(0)]], $P->[0]) if $n == 0;
  my @q; $q[$n-1] = $P->[$n]; for (my $b = $n - 1; $b >= 1; $b--) { $q[$b-1] = padd($P->[$b], pmul($r, $q[$b])) }
  my $rem = padd($P->[0], pmul($r, $q[0])); return (btrim(\@q), $rem) }
sub bev { my ($P, $j, $y) = @_; my $s = Q(0); my $yy = Q(1); for my $b (0..$#$P) { $s += pev($P->[$b], $j) * $yy; $yy *= $y } return $s }
# ---------- scalar series in t, truncated at order NT ----------
our $NT = 40;
sub szero { return [ map { Q(0) } 0..$NT ] }
sub smul { my ($a, $b) = @_; my $r = szero(); for my $i (0..$NT) { next if $a->[$i] == 0; for my $k (0..$NT-$i) { next if $b->[$k] == 0; $r->[$i+$k] += $a->[$i] * $b->[$k] } } return $r }
sub sadd { my ($a, $b, $c) = @_; return [ map { defined $c ? $a->[$_] + $b->[$_] * $c : $a->[$_] + $b->[$_] } 0..$NT ] }
sub ssc { my ($a, $c) = @_; return [ map { $_ * $c } @$a ] }
sub sexp { my $a = shift; die "sexp: nonzero constant term" if $a->[0] != 0; my $r = szero(); $r->[0] = Q(1);
  for my $n (1..$NT) { my $s = Q(0); for my $m (1..$n) { next if $a->[$m] == 0; $s += $a->[$m] * $r->[$n-$m] * $m } $r->[$n] = $s / $n } return $r }
sub sder { my $a = shift; my $r = szero(); for my $n (0..$NT-1) { $r->[$n] = $a->[$n+1] * ($n + 1) } return $r }
sub stdt { my $a = shift; return [ map { $a->[$_] * $_ } 0..$NT ] }   # t d/dt
# log S, S = (t/2)/sinh(t/2):  [t^{2m}] = -B_{2m}/(2m (2m)!)
sub logS { my $s = szero(); for my $m (1..int($NT/2)) { $s->[2*$m] = -bern(2*$m) / (Q(2*$m) * fact(2*$m)) } return $s }
# log H = log S - t (log S)'
sub logH { my $f = logS(); return sadd($f, stdt($f), Q(-1)) }
# polynomial in j of [t^n] (A(t) exp(j L(t))) = Sum_k j^k/k! [t^n] A L^k   (L = O(t^2))
sub coeff_expj { my ($A, $L, $n) = @_; my @p; my $Lk = szero(); $Lk->[0] = Q(1); my $k = 0;
  while (1) { my $prod = smul($A, $Lk); push @p, $prod->[$n] / fact($k); $k++; last if 2*$k > $n; $Lk = smul($Lk, $L) } return ptrim(\@p) }
1;
