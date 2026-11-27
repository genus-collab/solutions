#!/usr/bin/perl
# ref_flow.pl -- referee 2: independent checks of the ingredients of the T1 proof (s.4-s.7), mod p.
#   branch J tau        : d/dl [u^-j]_+ (1+tau) = j Pi_j, j <= J                                   (Lemma 4.1)
#   res  N K tau        : Res_mu u_t(mu)^n = (-1)^(n-1) t^(n-1) e^(-n tau^2 t) Pi_n H_n(1-e^-t), n <= N, t-order <= K
#                         (u_t = u e^(-K_t), K_t = sum_k kappa_k t^k ((1+tau) mu^k - tau (mu-1)^k))      (s.7(c))
#   curve K tau         : on the leaf curve, tau log((phi-1)/(l-1)) - (1+tau) log(phi/l) - K_t(phi) = -log C,
#                         C = W e^(tau^2 t)/H_1, independent of l                                     (Theorem 7.1)
#   path K tau seed dens: along a random path P_j(t) = t^(j-1)(random), F_G(P(t)) - F_G(0) from the KP engine (rk.pm,
#                         symbolic in the one-leaf values) versus Res_mu Psi(mu) [log W + log u(mu) - log u(psi_t(mu))]
#                         (Theorem 6.1, flow theorem)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rk;
my ($mode, @a) = @ARGV;
sub frac { my ($n, $d) = split m{/}, shift; $d //= 1; return mfrac($n, $d); }
sub Pis { my ($tau, $n) = @_; my @Pi = (1, 1); for my $j (2..$n) { my $v = 1; for my $b (1..$j-1) { $v = mm($v, md($j*$tau + $b)); } $Pi[$j] = mm($v, ifct($j-1)); } return \@Pi; }
# [u^-j]_+ coefficients: u^-j = l^j (1-1/l)^(-j tau) = sum_i binom(-j tau, i) (-1)^i l^(j-i)
sub uplus { my ($tau, $j) = @_; my @c = (0) x ($j+1); for my $i (0..$j) { $c[$j-$i] = mm(gbin(md(-$j*$tau), $i), ($i % 2) ? $P - 1 : 1); } return \@c; }
# Res_mu u^n mu^d = [mu^-1] mu^(d-n) (1-1/mu)^(n tau) = (-1)^(d-n+1) binom(n tau, d-n+1)
sub resu { my ($tau, $n, $d) = @_; my $e = $d - $n + 1; return 0 if $e < 0; return mm(gbin(mm($n, $tau), $e), ($e % 2) ? $P - 1 : 1); }
# ---- bivariate helpers: t-series (0..K) of polynomials (arrays) ----
sub bz { my $K = shift; return [ map { [0] } 0..$K ]; }
sub badd { my ($x, $y, $c, $K) = @_; $c //= 1; my @r; for my $k (0..$K) { my @p = @{$x->[$k]}; my $q = $y->[$k]; for my $d (0..$#$q) { $p[$d] = md(($p[$d] // 0) + mm($c, $q->[$d])); } $r[$k] = \@p; } return \@r; }
sub bmul { my ($x, $y, $K) = @_; my @r = map { [0] } 0..$K;
  for my $k1 (0..$K) { my $p = $x->[$k1]; next unless grep { $_ } @$p; for my $k2 (0..$K-$k1) { my $q = $y->[$k2]; next unless grep { $_ } @$q; my $o = $r[$k1+$k2];
      for my $i (0..$#$p) { next unless $p->[$i]; for my $j (0..$#$q) { next unless $q->[$j]; $o->[$i+$j] = md(($o->[$i+$j] // 0) + $p->[$i]*$q->[$j]); } } } }
  for my $o (@r) { $_ //= 0 for @$o; } return \@r; }
sub bscal { my ($x, $c) = @_; return [ map { [ map { mm($c, $_) } @$_ ] } @$x ]; }
# log(1 + X) for X with X[0] = 0
sub blog1 { my ($X, $K) = @_; my $r = bz($K); my $pw = [ [1], map { [0] } 1..$K ]; for my $n (1..$K) { $pw = bmul($pw, $X, $K); $r = badd($r, $pw, mm(($n % 2) ? 1 : $P - 1, mi($n)), $K); } return $r; }
# divide each polynomial by (l - r) exactly (synthetic division), dying if the remainder is non-zero
sub bdivlin { my ($X, $r) = @_; my @out; for my $p (@$X) { my @c = @$p; my $n = $#c; if ($n < 1) { die "bdivlin: remainder\n" if $c[0]; push @out, [0]; next; }
    my @q = (0) x $n; my $acc = 0; for (my $d = $n; $d >= 1; $d--) { $acc = md($c[$d] + mm($acc, $r)); $q[$d-1] = $acc; } my $rem = md($c[0] + mm($acc, $r)); die "bdivlin: remainder $rem\n" if $rem; push @out, \@q; } return \@out; }
# compose: f(l) = sum_k t^k f_k(l) evaluated at l = g (a t-series of polynomials in mu with g_0 = mu)
sub bcompose { my ($f, $g, $K) = @_; my $maxd = 0; for my $p (@$f) { $maxd = $#$p if $#$p > $maxd; }
  my @pw = ([ [1], map { [0] } 1..$K ]); for my $d (1..$maxd) { $pw[$d] = bmul($pw[$d-1], $g, $K); }
  my $r = bz($K); for my $k (0..$K) { my $p = $f->[$k]; for my $d (0..$#$p) { next unless $p->[$d]; my $sh = [ map { $_ < $k ? [0] : $pw[$d][$_ - $k] } 0..$K ]; $r = badd($r, $sh, $p->[$d], $K); } } return $r; }
# phi_t(l) = W l - [Pcal]_+(l) + [Pcal]_+(0), Pcal = sum_j Phat_j u^-j; given W (t-series) and Phat_j (t-series)
sub phi_of { my ($tau, $W, $Phat, $K) = @_; my $phi = bz($K); for my $k (0..$K) { $phi->[$k] = [0, $W->[$k]]; }
  for my $j (keys %$Phat) { my $up = uplus($tau, $j); my $c0 = $up->[0]; for my $k (0..$K) { my $s = $Phat->{$j}[$k]; next unless $s;
      my @p = @{$phi->[$k]}; for my $d (1..$j) { $p[$d] = md(($p[$d] // 0) - mm($s, $up->[$d])); } $phi->[$k] = [ map { $_ // 0 } @p ]; } }
  return $phi; }
sub inverse { my ($phi, $K) = @_; my $delta = [ map { [ @$_ ] } @$phi ]; $delta->[0][1] = md($delta->[0][1] - 1); die "phi_0 != l\n" if grep { $_ } @{$delta->[0]};
  my $mu = bz($K); $mu->[0] = [0, 1]; my $psi = $mu;
  for my $it (1..$K) { my $dp = bcompose($delta, $psi, $K); $psi = badd($mu, $dp, $P - 1, $K); } return $psi; }
if ($mode eq 'branch') { my ($J, $ts) = @a; my $tau = frac($ts); my $Pi = Pis($tau, $J); my $bad = 0; my $l0 = md(1 + $tau);
  for my $j (1..$J) { my $up = uplus($tau, $j); my $v = 0; my $pw = 1; for my $d (1..$j) { $v = md($v + mm(mm($d, $up->[$d]), $pw)); $pw = mm($pw, $l0); } $bad++ if $v != mm($j, $Pi->[$j]); }
  print "branch tau=$ts: j <= $J, $bad mismatches\n"; exit; }
if ($mode eq 'res') { my ($N, $K, $ts) = @a; my $tau = frac($ts); my $Pi = Pis($tau, $N + 2); my $kap = kappa_list($tau, $K); my $bad = 0; my $cnt = 0;
  # K_t as t-series of polynomials in mu
  my $Kt = bz($K); for my $k (1..$K) { my @p = (0) x ($k+1); $p[$k] = md(1 + $tau); for my $i (0..$k) { $p[$i] = md($p[$i] - mm($tau, mm(binom_int($k, $i), (($k - $i) % 2) ? $P - 1 : 1))); } $Kt->[$k] = [ map { mm($kap->[$k], $_) } @p ]; }
  my $T = T_series($K);
  for my $n (1..$N) { my $E = bz($K); $E->[0] = [1]; my $pw = [ [1], map { [0] } 1..$K ]; my $X = bscal($Kt, md(-$n));
    for my $m (1..$K) { $pw = bmul($pw, $X, $K); $E = badd($E, $pw, ifct($m), $K); }      # e^{-n K_t}
    my @lhs = map { my $p = $E->[$_]; my $s = 0; for my $d (0..$#$p) { $s = md($s + mm($p->[$d], resu($tau, $n, $d))); } $s } 0..$K;
    my $H = ser_compose(hyp_series(md(-$n*$tau), md(-$n*$tau), $n, $K), $T, $K); my @ex = map { mm(ifct($_), mpow(md(-mm($n, mm($tau, $tau))), $_)) } 0..$K;
    my $rhs = ser_mul($H, \@ex, $K); my $c = mm(($n % 2) ? 1 : $P - 1, $Pi->[$n]);
    for my $k (0..$K) { my $want = $k >= $n - 1 ? mm($c, $rhs->[$k - $n + 1]) : 0; $cnt++; $bad++ if $want != $lhs[$k]; } }
  print "res tau=$ts: n <= $N, t-order <= $K: $cnt coefficients, $bad mismatches\n"; exit; }
if ($mode eq 'curve') { my ($K, $ts) = @a; my $tau = frac($ts); my $Pi = Pis($tau, $K + 3); my $J = $K + 2;
  my $Pj = leaf_series($tau, $Pi, $K, $J);
  # W = 1 + sum Pi_j P_j W^j
  my $W = [1, (0) x $K]; for my $it (1..$K+1) { my $s = [1, (0) x $K]; my $wp = [1, (0) x $K]; for my $j (1..$J) { $wp = ser_mul($wp, $W, $K); next if $j < 2; my $tj = ser_mul($Pj->{$j}, $wp, $K); $s->[$_] = md($s->[$_] + mm($Pi->[$j], $tj->[$_])) for 0..$K; } $W = $s; }
  my %Phat; my $wp = [1, (0) x $K]; for my $j (1..$J) { $wp = ser_mul($wp, $W, $K); next if $j < 2; $Phat{$j} = ser_mul($Pj->{$j}, $wp, $K); }
  my $phi = phi_of($tau, $W, \%Phat, $K);
  # check phi(0) = 0, phi(1) = 1
  my $bad01 = 0; for my $k (0..$K) { my $p = $phi->[$k]; my $v1 = 0; $v1 = md($v1 + $_) for @$p; $bad01++ if md($p->[0]) != 0 || $v1 != ($k == 0 ? 1 : 0); }
  my $lmu = bz($K); $lmu->[0] = [0, 1];
  my $A = blog1(bdivlin(badd($phi, $lmu, $P - 1, $K), 0), $K);   # log(phi/l): (phi - l)/l
  my $Bv = blog1(bdivlin(badd($phi, $lmu, $P - 1, $K), 1), $K);  # log((phi-1)/(l-1)): (phi - l)/(l - 1)
  my $kap = kappa_list($tau, $K); my $Kt = bz($K);
  for my $k (1..$K) { my @p = (0) x ($k+1); $p[$k] = md(1 + $tau); for my $i (0..$k) { $p[$i] = md($p[$i] - mm($tau, mm(binom_int($k, $i), (($k - $i) % 2) ? $P - 1 : 1))); } $Kt->[$k] = [ map { mm($kap->[$k], $_) } @p ]; }
  my $Kphi = bcompose($Kt, $phi, $K);
  my $Z = badd(badd(bscal($Bv, $tau), $A, md(-(1 + $tau)), $K), $Kphi, $P - 1, $K);
  # predicted: -log C = -log W - tau^2 t + log H_1
  my $T = T_series($K); my $H1 = ser_compose(hyp_series(md(-$tau), md(-$tau), 1, $K), $T, $K);
  my $lW = ser_log($W, $K); my $lH = ser_log($H1, $K); my @pred = map { md($lH->[$_] - $lW->[$_]) } 0..$K; $pred[1] = md($pred[1] - mm($tau, $tau));
  my $bad = 0; for my $k (0..$K) { my $p = $Z->[$k]; for my $d (1..$#$p) { $bad++ if $p->[$d]; } $bad++ if md($p->[0] // 0) != $pred[$k]; }
  print "curve tau=$ts K=$K: phi(0)=0, phi(1)=1 failures $bad01; curve identity (l-independence and value -log C) mismatches $bad\n"; exit; }
if ($mode eq 'path') { my ($K, $ts, $seed, $dens) = @a; $dens //= 1; srand($seed); my $tau = frac($ts); my $Pi = Pis($tau, 2*$K + 4); my $J = $K + 1;
  my %Pj; for my $j (2..$J) { next if rand() > $dens && $j > 2; my @s = (0) x ($K+1); for my $o ($j-1..$K) { $s[$o] = int(rand(2001)) - 1000; $s[$o] = md($s[$o]); } $Pj{$j} = \@s; }
  for my $j (2..$J) { $Pj{$j} //= [ (0) x ($K+1) ]; }
  # engine: symbolic one-leaf unknowns y_1..y_K
  my $m = $K; my @yv; for my $L (1..$K) { my @v = (0) x $m; $v[$L-1] = 1; $yv[$L] = \@v; }
  my $S = kp_setup($tau, $K, \@yv, 'big');
  my %PS; my @Ts = sort { scalar(@$a) <=> scalar(@$b) } @{$S->{Ts}};
  for my $T (@Ts) { my @t = @$T; my $s = pop @t; my $mult = scalar grep { $_ == $s } @$T; my $base = @t ? $PS{keyof(@t)} : [1, (0) x $K];
    my $pr = ser_mul($base, $Pj{$s}, $K); my $im = mi($mult); $PS{keyof(@$T)} = [ map { mm($_, $im) } @$pr ]; }
  my @Phi = map { [ (0) x $m ] } 0..$K;
  for my $T (@Ts) { my $k = keyof(@$T); my $ps = $PS{$k}; my $y = $S->{Y}{$k}; for my $o (0..$K) { my $c = $ps->[$o]; next unless $c; for my $i (0..$m-1) { $Phi[$o][$i] = ($Phi[$o][$i] + $c*$y->[$i]) % $P if $y->[$i]; } } }
  # closed formula
  my $W = [1, (0) x $K]; for my $it (1..$K+1) { my $s = [1, (0) x $K]; my $wp = [1, (0) x $K]; for my $j (1..$J) { $wp = ser_mul($wp, $W, $K); next if $j < 2; my $tj = ser_mul($Pj{$j}, $wp, $K); $s->[$_] = md($s->[$_] + mm($Pi->[$j], $tj->[$_])) for 0..$K; } $W = $s; }
  my %Phat; my $wp = [1, (0) x $K]; for my $j (1..$J) { $wp = ser_mul($wp, $W, $K); next if $j < 2; $Phat{$j} = ser_mul($Pj{$j}, $wp, $K); }
  my $phi = phi_of($tau, $W, \%Phat, $K); my $psi = inverse($phi, $K);
  my $lmu = bz($K); $lmu->[0] = [0, 1];
  my $A = blog1(bdivlin(badd($psi, $lmu, $P - 1, $K), 0), $K);   # log(psi/mu)
  my $Bv = blog1(bdivlin(badd($psi, $lmu, $P - 1, $K), 1), $K);  # log((psi-1)/(mu-1))
  # Kflow = log W + log u(mu) - log u(psi) = log W - tau log((psi-1)/(mu-1)) + (1+tau) log(psi/mu)
  my $Kf = badd(bscal($A, md(1 + $tau)), $Bv, md(-$tau), $K); my $lW = ser_log($W, $K); for my $k (0..$K) { $Kf->[$k][0] = md(($Kf->[$k][0] // 0) + $lW->[$k]); }
  my ($bad, $cnt, $nz) = (0, 0, 0);
  for my $k (0..$K) { my @pred = (0) x $m; my $p = $Kf->[$k];
    for my $d (0..$#$p) { next unless $p->[$d]; for my $c (2..$d+1) { my $r = resu($tau, $c, $d); next unless $r; die "level $c-1 > K\n" if $c - 1 > $K; $pred[$c-2] = md($pred[$c-2] + mm($p->[$d], $r)); } }
    for my $i (0..$m-1) { $cnt++; $bad++ if $pred[$i] != $Phi[$k][$i]; $nz++ if $Phi[$k][$i]; } }
  print "path tau=$ts K=$K seed=$seed dens=$dens: ", scalar(grep { grep { $_ } @{$Pj{$_}} } keys %Pj), " non-zero P_j; $cnt coefficients ($nz non-zero), $bad mismatches\n"; exit; }
die "mode?\n";
