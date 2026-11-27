# claimE_proof_checks.pl : exact checks of the three identities used in the proof of Claim E (report s.4, Theorem K1):
# (V1) the formal Bernoulli identity  X(w/(1+w)) - X(w) = -log(1+w) + w/2 + w/(2(1+w)),  X(w) = Sum_{m>=1} chi_m w^{2m}, chi_m = -B_{2m}/(2m);
# (V2) the generating function  Sum_{M>=1} C_ext(M;a,n) w^{2M+2} = X(w)[w^2 (a-1)(a+2n) + 2 - (1+w)^a - (1-w)^a]
#        + (2a-1) w^2 - (w/2)[(1+w)^a - (1-w)^a + (1+w)^{a-1} - (1-w)^{a-1}],  with
#      C_ext(M;a,n) = chi_M [2M(2a-2M-1) + 2n(a-1)] - 2 d/da binom(a,2M+2) - 2 Sum_{m=1}^{M-1} binom(a-2m, 2M+2-2m) chi_m
#      (checked coefficientwise, a = 1..12, n = 0..4, M <= 14);
# (V3) (R_T):  (N-1) Gamma_{N-1}(j) = Sum_{n=0}^{N-2} (j-2n) Gamma_n(j)/(2N-2n)!  as polynomials in j, Gamma_n = [t^{2n}] H^j
#      (equivalently t (log H)' = 1 - S^2), N <= 16.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use PJ;
sub Q { Math::BigRat->new(shift) }
my $W = 30;   # order in w
my @F = (Q(1)); sub fact { my $n = shift; while (@F <= $n) { push @F, $F[-1] * scalar(@F) } return $F[$n] }
sub binq { my ($x, $n) = @_; return Q(0) if $n < 0; my $p = Q(1); for my $i (0..$n-1) { $p *= ($x - $i) } return $p / fact($n) }
sub dbinq { my ($x, $n) = @_; my $s = Q(0); for my $i (0..$n-1) { my $p = Q(1); for my $l (0..$n-1) { next if $l == $i; $p *= ($x - $l) } $s += $p } return $s / fact($n) }
sub chi { my $m = shift; return -PJ::bern(2*$m) / (2*$m) }
# series in w as arrays 0..W
sub zero { [ map { Q(0) } 0..$W ] }
sub add { my ($a, $b, $c) = @_; $c //= Q(1); [ map { $a->[$_] + $b->[$_] * $c } 0..$W ] }
sub mul { my ($a, $b) = @_; my $r = zero(); for my $i (0..$W) { next if $a->[$i] == 0; for my $k (0..$W-$i) { $r->[$i+$k] += $a->[$i] * $b->[$k] } } $r }
sub compose { my ($A, $B) = @_; # A(B(w)), B(0) = 0
  my $r = zero(); my $p = zero(); $p->[0] = Q(1); for my $k (0..$W) { $r = add($r, $p, $A->[$k]) if $A->[$k] != 0; $p = mul($p, $B); } $r }
my $X = zero(); for my $m (1..int($W/2)) { $X->[2*$m] = chi($m) }
# (V1)
my $u = zero(); for my $k (1..$W) { $u->[$k] = Q((-1) ** ($k - 1)) }   # w/(1+w) = Sum (-1)^{k-1} w^k
my $lhs = add(compose($X, $u), $X, Q(-1));
my $rhs = zero(); for my $k (1..$W) { $rhs->[$k] -= Q((-1) ** ($k - 1)) / $k; }   # -log(1+w)
$rhs->[1] += Q(1)/2; for my $k (1..$W) { $rhs->[$k] += Q((-1) ** ($k - 1)) / 2; }   # w/2 + (w/2)/(1+w)
my $bad1 = grep { $lhs->[$_] != $rhs->[$_] } 0..$W;
print "(V1) formal Bernoulli identity to order w^$W: ", ($bad1 ? "$bad1 FAILURES" : "holds"), "\n";
# (V2)
my $bad2 = 0; my $n2 = 0;
for my $a (1..12) { for my $n (0..4) {
    my $direct = zero();
    for my $M (1..14) { last if 2*$M + 2 > $W; my $c = chi($M) * (2*$M*(2*$a - 2*$M - 1) + 2*$n*($a - 1)) - 2 * dbinq(Q($a), 2*$M + 2);
      for my $m (1..$M-1) { $c -= 2 * binq(Q($a - 2*$m), 2*$M + 2 - 2*$m) * chi($m) } $direct->[2*$M + 2] = $c; }
    my $p1 = zero(); my $p2 = zero(); for my $k (0..$W) { $p1->[$k] = binq(Q($a), $k); $p2->[$k] = binq(Q($a), $k) * Q((-1) ** $k); }   # (1+w)^a, (1-w)^a
    my $q1 = zero(); my $q2 = zero(); for my $k (0..$W) { $q1->[$k] = binq(Q($a - 1), $k); $q2->[$k] = binq(Q($a - 1), $k) * Q((-1) ** $k); }
    my $br = add(add(zero(), $p1, Q(-1)), $p2, Q(-1)); $br->[0] += 2; $br->[2] += ($a - 1) * ($a + 2*$n);
    my $cf = mul($X, $br); $cf->[2] += 2*$a - 1;
    my $d = add(add(add($p1, $p2, Q(-1)), $q1), $q2, Q(-1)); for my $k (0..$W-1) { $cf->[$k+1] -= $d->[$k] / 2 }
    for my $k (0..$W) { next if $k <= 2*14 + 2 && $k >= 0 && 0; }
    for my $M (1..14) { my $k = 2*$M + 2; last if $k > $W; $n2++; $bad2++ if $cf->[$k] != $direct->[$k]; } } }
print "(V2) generating function of C_ext: $n2 coefficients compared, $bad2 failures\n";
# (V3)
$PJ::N = 34; my $f = PJ::logS(); my $logH = PJ::sadd($f, PJ::stder($f), Q(-1)); my $G = PJ::sexp(PJ::sscalp($logH, PJ::pj()));
my $bad3 = 0;
for my $N (2..16) { my $l = PJ::pscale($G->[2*$N - 2], Q($N - 1)); my $r = PJ::pz();
  for my $n (0..$N-2) { $r = PJ::padd($r, PJ::pscale(PJ::pmul([Q(-2*$n), Q(1)], $G->[2*$n]), Q(1) / fact(2*$N - 2*$n))); }
  $bad3++ unless PJ::iszero(PJ::psub($l, $r)); }
print "(V3) (R_T) as polynomial identities in j for N = 2..16: ", ($bad3 ? "$bad3 FAILURES" : "holds"), "\n";
