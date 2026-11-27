#!/usr/bin/perl
# moments_check.pl tau kmax : Laplace moments behind the dictionary P2 (Eynard/DOSS, as quoted in ABDKS arXiv:2312.16950 s.3.5, (49)-(50)).
# On C_tau (X = tau log v - e2 log(1+v), y = -v), the leaves xi_0 = 1/(v - tau), xi_k = D^k xi_0, D = d/dX = (v(1+v)/(tau - v)) d/dv.
# q_k is proportional to Res_{v=tau} y dxi_k = Res_{v=tau} xi_k dv (integration by parts); the claim is Res_{v=tau} xi_k dv = 1 for all k,
# i.e. q_k/q_0 = (+-1)^k, which gives exp(Sum s_j kappa_j) = exp(-Sum kappa_j/j) = K.
# Method: write xi_k = N_k(v)/(v - tau)^{2k+1} with N_k a polynomial; compute residues exactly via Laurent expansion at v = tau.
use strict; use warnings; use Math::BigRat;
sub Q { Math::BigRat->new(@_) }
my ($tau, $kmax) = @ARGV; $tau = Q($tau // 2); $kmax //= 6;
# polynomials in u = v - tau as arrays (index = power)
sub padd { my ($a, $b) = @_; my $n = $#$a > $#$b ? $#$a : $#$b; return [ map { ($a->[$_] // Q(0)) + ($b->[$_] // Q(0)) } 0..$n ]; }
sub pmul { my ($a, $b) = @_; my @c = map { Q(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { for my $j (0..$#$b) { $c[$i+$j] += $a->[$i] * $b->[$j]; } } return \@c; }
sub pscal { my ($s, $a) = @_; return [ map { $_ * $s } @$a ]; }
sub pder { my ($a) = @_; return [Q(0)] if $#$a < 1; return [ map { $a->[$_] * $_ } 1..$#$a ]; }
# v = tau + u;  v(1+v) = (tau + u)(tau + 1 + u)
my $V1 = pmul([ $tau, Q(1) ], [ $tau + 1, Q(1) ]);
# xi_k = N_k / u^{m_k}; D = -(v(1+v)/u) d/du ; D(N/u^m) = -(V1/u) (N'/u^m - m N/u^{m+1}) = (-V1 N' u + m V1 N)/u^{m+2}
my ($N, $m) = ([Q(1)], 1);
for my $k (0..$kmax) {
  # residue at u = 0 of N/u^m: coefficient of u^{m-1} in N
  my $res = $N->[$m-1] // Q(0);
  printf "tau=%s k=%d: Res_{v=tau} xi_k dv = %s\n", $tau, $k, $res;
  my $t1 = pscal(Q(-1), pmul(pmul($V1, pder($N)), [Q(0), Q(1)]));
  my $t2 = pscal(Q($m), pmul($V1, $N));
  $N = padd($t1, $t2); $m += 2; }
