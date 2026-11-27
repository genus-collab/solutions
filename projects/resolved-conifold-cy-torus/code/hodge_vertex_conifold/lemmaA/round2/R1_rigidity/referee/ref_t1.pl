#!/usr/bin/perl
# ref_t1.pl -- referee's independent test of T1 (closed form of the polar-order relations).
#   relation_k := [t^k] Sum_T prod_{t in T} P_t(t) Y(T)/|Aut T|,  P_t = (-1)^t (Pi_t/t) t^(t-1) H_t(T)/H_1(T)^t,
#   H_j = 2F1(-j tau, -j tau; j; T), T = 1 - e^-t   (leaf curve: Theorem 11.1 + string elimination),
# with Y(T) propagated from arbitrary one-leaf values by refkp::kp_propagate (my own anti-diagonal form of the
# linearised dispersionless Hirota identity at every point of the slice), compared with
#   T1: relation_k = kappa_k Sum_L y_L (-1)^(k-L) [ (1+tau) binom((L+1)tau, k-L) - tau binom((L+1)tau + k, k-L) ],
#   kappa_k = [x^k] (tau x + log((e^x-1)/x))   (computed from the series, not from Bernoulli numbers).
# Usage: perl ref_t1.pl mode(mod|rat) K tau1 tau2 ...   [env PRIME=...]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use refkp;
my ($mode, $K, @taus) = @ARGV; set_mode($mode); set_prime($ENV{PRIME}) if $ENV{PRIME};
my ($tot, $bad, $chk, $chkbad) = (0, 0, 0, 0);
for my $ts (@taus) {
  my ($a, $b) = split m{/}, $ts; $b //= 1; my $tau = FR($a, $b);
  my $t0 = time; my $S = kp_propagate($tau, $K); $chk += $S->{nchk}; $chkbad += $S->{nbad};
  my $Pi = $S->{Pi};
  # leaf curve series
  my $N = $K; my $Tt = Tser($N);
  my $h1 = ser_comp(hyp2f1(NG($tau), NG($tau), ONE(), $N), $Tt, $N);
  my %rho;
  for my $j (2..$K+1) { my $hj = ser_comp(hyp2f1(NG(ML(NN($j), $tau)), NG(ML(NN($j), $tau)), NN($j), $N), $Tt, $N);
    my $pw = [ ONE(), map { ZERO() } 1..$N ]; $pw = ser_mul($pw, $h1, $N) for 1..$j; $rho{$j} = ser_mul($hj, ser_inv($pw, $N), $N); }
  my @rel = map { vzero($K) } 0..$K;
  for my $D (1..$K) { for my $T (parts_ge2($D, $D + 1)) {
      my $top = ONE(); my $pr = [ ONE(), map { ZERO() } 1..$K ];
      for my $t (@$T) { $top = ML($top, ML(NN(($t % 2) ? -1 : 1), ML($Pi->($t), IV(NN($t))))); $pr = ser_mul($pr, $rho{$t}, $K - $D); }
      $top = ML($top, IV(aut_f(@$T)));
      my $v = $S->{Y}{mkey(@$T)}; die "no Y(@$T)" unless $v;
      for my $k ($D..$K) { my $c = ML($top, $pr->[$k - $D]); vaxpy($rel[$k], $c, $v); } } }
  # kappa_k from the series log(sum x^n/(n+1)!)
  my @ex = map { IV(fact_f($_ + 1)) } 0..$K; my $lg = ser_log(\@ex, $K); my @kap = @$lg; $kap[1] = AD($kap[1], $tau);
  my $one = ONE();
  for my $k (1..$K) { for my $L (1..$K) { my $pred = ZERO();
      if ($L <= $k) { my $Nn = $L + 1; my $d = $k - $L;
        my $g1 = binom_gen(ML(NN($Nn), $tau), $d); my $g2 = binom_gen(AD(ML(NN($Nn), $tau), NN($k)), $d);
        $pred = ML($kap[$k], ML(NN(($d % 2) ? -1 : 1), SB(ML(AD($one, $tau), $g1), ML($tau, $g2)))); }
      $tot++; unless (EQ($pred, $rel[$k][$L-1])) { $bad++; print "MISMATCH tau=$ts k=$k L=$L: got $rel[$k][$L-1] pred $pred\n" if $bad <= 10; } } }
  printf "tau=%s K=%d: propagation %d consistency checks (%d failed), time %ds; kappa_2..4 = %s %s %s\n", $ts, $K, $S->{nchk}, $S->{nbad}, time - $t0,
    map { "$kap[$_]" } 2..4;
  if ($ENV{SHOW}) { for my $k (1..($K < 4 ? $K : 4)) { print "  relation_$k: ", join(' ', map { "$_" } @{$rel[$k]}), "\n"; } }
}
print "ref_t1 ($mode): $tot coefficients compared, $bad mismatches; propagation checks $chk, failed $chkbad\n";
