#!/usr/bin/perl
# ref_res.pl -- referee's check of the residue identity of s.7(c):
#   Res_mu u_t(mu)^n dmu = (-1)^(n-1) t^(n-1) e^(-n tau^2 t) Pi_n 2F1(-n tau, -n tau; n; 1 - e^-t),
#   u_t^n = u^n exp(-n K_t),  K_t(mu) = Sum_k kappa_k t^k ((1+tau) mu^k - tau (mu-1)^k),  Res u^n mu^i = (-1)^(i-n+1) binom(n tau, i-n+1).
# kappa_k from the series of log((e^x-1)/x) (+ tau x).  Usage: perl ref_res.pl M nmax tau...   (mod p)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use refkp;
my ($M, $nmax, @taus) = @ARGV; set_mode('mod'); my $one = ONE(); my ($tot, $bad) = (0, 0);
for my $ts (@taus) { my ($a, $b) = split m{/}, $ts; $b //= 1; my $tau = FR($a, $b);
  my @ex = map { IV(fact_f($_ + 1)) } 0..$M; my $lg = ser_log(\@ex, $M); my @kap = @$lg; $kap[1] = AD($kap[1], $tau);
  # K[k] = polynomial in mu (array) for t^k
  my @K = ([ZERO()]); for my $k (1..$M) { my @p = map { ZERO() } 0..$k; $p[$k] = AD($p[$k], AD($one, $tau));
    for my $i (0..$k) { $p[$i] = SB($p[$i], ML($tau, ML(binom_gen(NN($k), $i), NN((($k-$i) % 2) ? -1 : 1)))); } $K[$k] = [ map { ML($kap[$k], $_) } @p ]; }
  for my $n (1..$nmax) {
    # E = exp(-n K) as t-series of polynomials: E_m = (1/m) Sum_{k=1}^m k (-n K_k) E_{m-k}
    my @E = ([$one]);
    for my $m (1..$M) { my @s = map { ZERO() } 0..$m; for my $k (1..$m) { my $A = $K[$k]; my $B = $E[$m-$k];
        for my $i (0..$#$A) { next if IZ($A->[$i]); for my $j (0..$#$B) { $s[$i+$j] = AD($s[$i+$j], ML(ML(NN(-$n*$k), $A->[$i]), $B->[$j])); } } }
      $E[$m] = [ map { ML($_, IV(NN($m))) } @s ]; }
    my @lhs = map { my $m = $_; my $s = ZERO(); for my $i (0..$#{$E[$m]}) { my $e = $i - $n + 1; next if $e < 0; $s = AD($s, ML($E[$m][$i], ML(NN(($e % 2) ? -1 : 1), binom_gen(ML(NN($n), $tau), $e)))); } $s } 0..$M;
    my $Tt = Tser($M); my $hn = ser_comp(hyp2f1(NG(ML(NN($n), $tau)), NG(ML(NN($n), $tau)), NN($n), $M), $Tt, $M);
    my $ee = ser_exp([ ZERO(), NG(ML(NN($n), ML($tau, $tau))), map { ZERO() } 2..$M ], $M);
    my $r = ser_mul($hn, $ee, $M); my $Pn = Pi_of($n, $tau); my $sg = NN((($n-1) % 2) ? -1 : 1);
    my @rhs = map { ZERO() } 0..$M; for my $m ($n-1..$M) { $rhs[$m] = ML($sg, ML($Pn, $r->[$m-$n+1])); }
    for my $m (0..$M) { $tot++; $bad++ unless EQ($lhs[$m], $rhs[$m]); } } }
print "ref_res: $tot t-coefficients compared, $bad mismatches\n";
