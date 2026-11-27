#!/usr/bin/perl
# ref_flow.pl -- referee's independent tests of the flow theorem (Theorem 6.1) and of the curve identity (Theorem 7.1).
# (A) Random path P_j(t) = O(t^(j-1)) in the slice.  Direct side: Phi(t) = Sum_T Y(T) P(t)^T/|Aut T| with Y(T) from
#     refkp::kp_propagate (one-leaf values = free unknowns).  Closed side (route's s.6):
#        Phi(t) = Res_mu Psi(mu) [ log W + log u(mu) - log u(psi_t(mu)) ],  psi_t = phi_t^{-1},
#        phi_t(l) = W l - [Pcal]_+(l) + [Pcal]_+(0),  Pcal = Sum_j P_j W^j u^{-j},  W = 1 + Sum_j Pi_j P_j W^j,
#        [u^{-j}]_+(l) = Sum_{k=0}^{j} binom(j tau + k - 1, k) l^(j-k),  Res u^c mu^n = (-1)^(n-c+1) binom(c tau, n-c+1).
#     Also checks phi_t(0) = 0, phi_t(1) = 1.
# (B) Leaf curve: Kflow_t(mu) - K_t(mu) must be independent of mu (equal to log(alpha/C) = log H_1 - tau^2 t),
#     K_t(mu) = Sum_k kappa_k t^k ((1+tau) mu^k - tau (mu-1)^k); and W = alpha := zeta H_1/t equals the root of the genus-0 equation.
# Usage: perl ref_flow.pl K tau [seed]   (mod p)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use refkp;
my ($K, $ts, $seed) = @ARGV; $seed //= 7; set_mode('mod'); set_prime($ENV{PRIME}) if $ENV{PRIME};
my ($ta, $tb) = split m{/}, $ts; $tb //= 1; my $tau = FR($ta, $tb); my $one = ONE();
my $DEG = $K + 3;
# bivariate: B->[m][i] coefficient of t^m mu^i, m = 0..K, i = 0..DEG
sub bz { return [ map { [ map { ZERO() } 0..$DEG ] } 0..$K ]; }
sub bmul { my ($A, $B) = @_; my $C = bz();
  for my $m1 (0..$K) { for my $i1 (0..$DEG) { my $a = $A->[$m1][$i1]; next if IZ($a);
      for my $m2 (0..$K-$m1) { for my $i2 (0..$DEG-$i1) { my $b = $B->[$m2][$i2]; next if IZ($b); $C->[$m1+$m2][$i1+$i2] = AD($C->[$m1+$m2][$i1+$i2], ML($a, $b)); } } } }
  return $C; }
sub badd { my ($A, $B, $c) = @_; $c //= $one; my $C = bz(); for my $m (0..$K) { for my $i (0..$DEG) { $C->[$m][$i] = AD($A->[$m][$i], ML($c, $B->[$m][$i])); } } return $C; }
sub bconst_t { my $s = shift; my $C = bz(); for my $m (0..$K) { $C->[$m][0] = $s->[$m] // ZERO(); } return $C; }  # t-series times mu^0
sub bmu { my $C = bz(); $C->[0][1] = $one; return $C; }
# compose: f(l) = Sum_i f_i(t) l^i (bivariate in (t,l)) evaluated at l = g(t,mu)
sub bcomp { my ($f, $g) = @_; my $res = bz(); my $pw = bz(); $pw->[0][0] = $one;
  for my $i (0..$DEG) { my $fi = bz(); my $nz = 0; for my $m (0..$K) { $fi->[$m][0] = $f->[$m][$i]; $nz ||= !IZ($f->[$m][$i]); }
    if ($nz) { $res = badd($res, bmul($fi, $pw)); }
    $pw = bmul($pw, $g); } return $res; }
# log(1+X), X = O(t)
sub blog1p { my $X = shift; my $res = bz(); my $pw = bz(); $pw->[0][0] = $one;
  for my $n (1..$K) { $pw = bmul($pw, $X); my $c = ML(NN(($n % 2) ? 1 : -1), IV(NN($n))); $res = badd($res, $pw, $c); } return $res; }
# divide by mu (requires zero constant term) and by (mu - 1) (requires value 0 at mu = 1)
sub bdivmu { my $A = shift; my $C = bz(); my $ok = 1; for my $m (0..$K) { $ok &&= IZ($A->[$m][0]); for my $i (1..$DEG) { $C->[$m][$i-1] = $A->[$m][$i]; } } return ($C, $ok); }
sub bdivmum1 { my $A = shift; my $C = bz(); my $ok = 1;
  for my $m (0..$K) { # synthetic division by (mu - 1): A = (mu-1) Q + r
    my @a = @{$A->[$m]}; my @q = map { ZERO() } 0..$DEG; my $carry = ZERO();
    for (my $i = $DEG; $i >= 1; $i--) { $carry = AD($a[$i], $carry); $q[$i-1] = $carry; }
    my $r = AD($a[0], $carry); $ok &&= IZ($r); $C->[$m] = \@q; } return ($C, $ok); }

my $S = kp_propagate($tau, $K); my $Pi = $S->{Pi};
print "propagation: $S->{nchk} checks, $S->{nbad} failed\n";
my $NJ = $K + 1;
sub run_path {
  my ($Pser, $label, $leaf) = @_;    # Pser->{j} = t-series
  # W
  my @W = ($one, map { ZERO() } 1..$K);
  for my $it (0..$K+1) { my @nw = ($one, map { ZERO() } 1..$K); my $wp = [@W];
    for my $j (2..$NJ) { $wp = ser_mul($wp, \@W, $K); my $q = ser_mul([ map { ML($Pi->($j), $_) } @{$Pser->{$j}} ], $wp, $K); $nw[$_] = AD($nw[$_], $q->[$_]) for 0..$K; }
    @W = @nw; }
  # phi(t,l)
  my $phi = bz(); for my $m (0..$K) { $phi->[$m][1] = $W[$m]; }
  my @Wp = ([ $one, map { ZERO() } 1..$K ]); for my $j (1..$NJ) { $Wp[$j] = ser_mul($Wp[$j-1], \@W, $K); }
  for my $j (2..$NJ) { my $c = ser_mul($Pser->{$j}, $Wp[$j], $K);
    for my $k (0..$j-1) { my $bc = binom_gen(AD(ML(NN($j), $tau), NN($k - 1)), $k); my $deg = $j - $k;   # l^(j-k), j-k >= 1 (the l^0 term cancels)
      for my $m (0..$K) { $phi->[$m][$deg] = SB($phi->[$m][$deg], ML($bc, $c->[$m])); } } }
  # checks phi(0) = 0, phi(1) = 1
  my ($ok0, $ok1) = (1, 1); for my $m (0..$K) { $ok0 &&= IZ($phi->[$m][0]); my $s = ZERO(); $s = AD($s, $_) for @{$phi->[$m]}; $ok1 &&= EQ($s, $m == 0 ? $one : ZERO()); }
  # inverse psi: psi = mu - delta(psi), delta = phi - l
  my $delta = bz(); for my $m (0..$K) { for my $i (0..$DEG) { $delta->[$m][$i] = $phi->[$m][$i]; } } $delta->[0][1] = SB($delta->[0][1], $one);
  my $psi = bmu(); for my $it (1..$K+1) { $psi = badd(bmu(), bcomp($delta, $psi), NG($one)); }
  # check phi(psi) = mu
  my $chk = badd(bcomp($phi, $psi), bmu(), NG($one)); my $inv_ok = 1; for my $m (0..$K) { for my $i (0..$DEG) { $inv_ok &&= IZ($chk->[$m][$i]); } }
  # Kflow = log W - tau log((psi-1)/(mu-1)) + (1+tau) log(psi/mu)
  my ($q1, $okA) = bdivmu($psi); my $pm1 = badd($psi, bz()); $pm1->[0][0] = SB($pm1->[0][0], $one); my ($q2, $okB) = bdivmum1($pm1);
  $q1->[0][0] = SB($q1->[0][0], $one); $q2->[0][0] = SB($q2->[0][0], $one);
  my $Kf = badd(badd(bz(), blog1p($q1), AD($one, $ENV{PERTURB} ? ML(NN(2),$tau) : $tau)), blog1p($q2), NG($tau));
  my $logW = ser_log(\@W, $K); for my $m (0..$K) { $Kf->[$m][0] = AD($Kf->[$m][0], $logW->[$m]); }
  if ($leaf) {
    # K_t(mu) = Sum_k kappa_k t^k ((1+tau) mu^k - tau (mu-1)^k)
    my @ex = map { IV(fact_f($_ + 1)) } 0..$K; my $lg = ser_log(\@ex, $K); my @kap = @$lg; $kap[1] = AD($kap[1], $tau);
    my $Kt = bz();
    for my $k (1..$K) { for my $i (0..$k) { my $bnm = binom_gen(NN($k), $i); my $c = ML(AD($one, $tau), $i == $k ? $one : ZERO());
        $c = SB($c, ML($tau, ML($bnm, NN((($k - $i) % 2) ? -1 : 1)))); $Kt->[$k][$i] = AD($Kt->[$k][$i], ML($kap[$k], $c)); } }
    my $dif = badd($Kf, $Kt, NG($one)); my $muindep = 1; for my $m (0..$K) { for my $i (1..$DEG) { $muindep &&= IZ($dif->[$m][$i]); } }
    # predicted constant log H_1 - tau^2 t
    my $Tt = Tser($K); my $h1 = ser_comp(hyp2f1(NG($tau), NG($tau), $one, $K), $Tt, $K); my $lh = ser_log($h1, $K);
    $lh->[1] = SB($lh->[1], ML($tau, $tau)); my $cst = 1; for my $m (0..$K) { $cst &&= EQ($dif->[$m][0], $lh->[$m]); }
    # W versus alpha = zeta H_1/t computed from Theorem 11.1's zeta (independent of the genus-0 equation):
    # x = -tau/(1+tau), c = 1/(1+tau); s* = 2/(c(1+sqrt(1-4xT/c^2))); dphi_T = s*/(1-T s*); zeta = T e^{dphi/c}  (series in T, then T = 1-e^-t)
    my $x = NG(ML($tau, IV(AD($one, $tau)))); my $c = IV(AD($one, $tau)); my $NT = $K + 1;
    my @u = (ZERO(), NG(ML(NN(4), ML($x, IV(ML($c, $c))))), map { ZERO() } 2..$NT); $u[0] = $one;
    my $lu = ser_log(\@u, $NT); my $sq = ser_exp([ map { ML($_, IV(NN(2))) } @$lu ], $NT);
    my @den = map { ML($c, $_) } @$sq; $den[0] = AD($den[0], $c); my $ss = ser_inv(\@den, $NT); $ss = [ map { ML(NN(2), $_) } @$ss ];
    my @Ts = (ZERO(), @{$ss}[0..$NT-1]); my @om = map { NG($_) } @Ts; $om[0] = $one; my $dd = ser_mul($ss, ser_inv(\@om, $NT), $NT);
    my @dphi = (ZERO()); for my $kk (1..$NT) { $dphi[$kk] = ML($dd->[$kk-1], IV(NN($kk))); }
    my $e = ser_exp([ map { ML($_, IV($c)) } @dphi ], $NT); my @zT = (ZERO(), @{$e}[0..$NT-1]);
    my $Tt1 = Tser($NT); my $zeta = ser_comp(\@zT, $Tt1, $NT);   # zeta(t), zeta = t + O(t^2)
    my @zt = map { $zeta->[$_+1] } 0..$K;                         # zeta/t
    my $alpha = ser_mul(\@zt, $h1, $K);
    my $wa = 1; for my $m (0..$K) { $wa &&= EQ($alpha->[$m], $W[$m]); }
    printf "  leaf curve: W (genus-0 root) = zeta H_1/t (Theorem 11.1 zeta): %s; Kflow - K_t independent of mu: %s; constant = log H_1 - tau^2 t: %s\n",
      $wa ? 'yes' : 'NO', $muindep ? 'yes' : 'NO', $cst ? 'yes' : 'NO';
  }
  # Res Psi Kflow: Psi = Sum_c y_{c-1} u^c, vector index (c-1)-1
  my @flow = map { vzero($K) } 0..$K;
  for my $m (1..$K) { for my $c (2..$K+1) { my $s = ZERO();
      for my $n ($c-1..$DEG) { my $e = $n - $c + 1; my $r = ML(NN(($e % 2) ? -1 : 1), binom_gen(ML(NN($c), $tau), $e)); $s = AD($s, ML($r, $Kf->[$m][$n])); }
      $flow[$m][$c-2] = AD($flow[$m][$c-2], $s); } }
  # direct side
  my @dir = map { vzero($K) } 0..$K;
  for my $D (1..$K) { for my $T (parts_ge2($D, $D + 1)) { my $pr = [ $one, map { ZERO() } 1..$K ]; $pr = ser_mul($pr, $Pser->{$_}, $K) for @$T;
      my $ia = IV(aut_f(@$T)); my $v = $S->{Y}{mkey(@$T)}; for my $m ($D..$K) { vaxpy($dir[$m], ML($ia, $pr->[$m]), $v); } } }
  my ($n, $bad) = (0, 0); for my $m (1..$K) { $n++; $bad++ unless veq($flow[$m], $dir[$m]); }
  printf "%s: phi(0)=0 %s, phi(1)=1 %s, inverse ok %s, psi/mu exact %s, (psi-1)/(mu-1) exact %s; Phi(t) orders 1..%d: %d mismatches\n",
    $label, map({ $_ ? 'yes' : 'NO' } ($ok0, $ok1, $inv_ok, $okA, $okB)), $K, $bad;
  return $bad;
}
# (A) random paths
srand($seed);
my $tb2 = 0;
for my $trial (1..2) { my %Ps; for my $j (2..$NJ) { my @s = map { ZERO() } 0..$K; for my $m ($j-1..$K) { $s[$m] = FR(int(rand(19)) - 9, int(rand(5)) + 1); } $Ps{$j} = \@s; }
  $tb2 += run_path(\%Ps, "random path $trial", 0); }
# perturbation control: the formula must fail if phi is altered -- emulate by perturbing W in a third path is not needed;
# instead compare with a deliberately wrong direct side (drop the 3-leaf terms) for sensitivity.
# (B) leaf curve
{ my %Ps; my $Tt = Tser($K); my $h1 = ser_comp(hyp2f1(NG($tau), NG($tau), $one, $K), $Tt, $K);
  for my $j (2..$NJ) { my $hj = ser_comp(hyp2f1(NG(ML(NN($j), $tau)), NG(ML(NN($j), $tau)), NN($j), $K), $Tt, $K);
    my $pw = [ $one, map { ZERO() } 1..$K ]; $pw = ser_mul($pw, $h1, $K) for 1..$j; my $rho = ser_mul($hj, ser_inv($pw, $K), $K);
    my $cc = ML(NN(($j % 2) ? -1 : 1), ML($Pi->($j), IV(NN($j)))); my @s = map { ZERO() } 0..$K; for my $m ($j-1..$K) { $s[$m] = ML($cc, $rho->[$m-$j+1]); } $Ps{$j} = \@s; }
  $tb2 += run_path(\%Ps, "leaf curve", 1); }
print "ref_flow K=$K tau=$ts: total mismatching orders $tb2\n";
