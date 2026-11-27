#!/usr/bin/perl
# t_ansatz.pl -- test of the "deformed curve" ansatz for the genus-G linearized solution along the leaf curve:
#   F^Y(P(t)) = [lambda^{-1}] Psi(lambda) log( z_t(lambda) u(lambda) ),
#   z_t(lambda) solves  alpha(t) Y_0(z/alpha(t)) - sum_{j>=2} P_j(t) z^j = lambda + tau     (z = lambda + O(1)),
# against the closed form (T1):  [t^k] (polynomial part of log(z_t u)) = kappa_k [ (1+tau) lambda^k - tau (lambda-1)^k ].
# Only the polynomial part in lambda matters (Psi = O(lambda^-2)).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M, $NEG) = @ARGV; $taus //= '2/7'; $M //= 4; $NEG //= 4; my $tau = R($taus); my $sig = $tau*(1+$tau);
my $OFF = $NEG + 2; my $TOPP = 2*$M + 3; my $LEN = $OFF + $TOPP;   # lambda powers from -OFF to TOPP
sub lz { [ map { R(0) } 0..$LEN ] }
sub lmul { my ($a, $b) = @_; my $c = lz(); for my $i (0..$LEN) { next if $a->[$i]->is_zero; for my $j (0..$LEN) { next if $b->[$j]->is_zero; my $k = $i + $j - $OFF; next if $k < 0 || $k > $LEN; $c->[$k] = $c->[$k] + $a->[$i]*$b->[$j]; } } return $c; }
sub ladd { my ($a, $b, $s) = @_; $s //= R(1); return [ map { $a->[$_] + $s*$b->[$_] } 0..$LEN ]; }
sub lmon { my ($p, $co) = @_; my $v = lz(); my $i = $p + $OFF; $v->[$i] = $co // R(1) if $i >= 0 && $i <= $LEN; return $v; }
sub bz { [ map { lz() } 0..$M ] }
sub bmul { my ($a, $b) = @_; my $c = bz(); for my $i (0..$M) { for my $j (0..$M-$i) { $c->[$i+$j] = ladd($c->[$i+$j], lmul($a->[$i], $b->[$j])); } } return $c; }
sub badd { my ($a, $b, $s) = @_; return [ map { ladd($a->[$_], $b->[$_], $s) } 0..$M ]; }
sub bconst { my $v = bz(); $v->[0] = lmon(0, $_[0]); return $v; }
# leaf curve P_j(t), alpha(t)
my $x = -$tau/(1+$tau); my $c = 1 + $x; my $N = $M + 2; my $Tt = T_of_t_series($N+1);
my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $N+1), $Tt, $N+1);
my %Pj; for my $j (2..$M+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $N+1), $Tt, $N+1);
  my $p = [ R(1), map { R(0) } 1..$N+1 ]; $p = series_mul($p, $h1, $N+1) for 1..$j; my $rho = series_mul($hj, series_inv($p, $N+1), $N+1);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$M ]; }
my @W = (R(1), map { R(0) } 1..$M);
for my $it (1..$M+1) { my @nw = (R(1), map { R(0) } 1..$M); my $wp = [@W];
  for my $j (2..$M+1) { $wp = series_mul($wp, \@W, $M); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp, $M); $nw[$_] = $nw[$_] + $q->[$_] for 0..$M; } @W = @nw; }
# z_0(lambda) = 1/u(lambda) = lambda (1 - 1/lambda)^{-tau}
my $z0 = lz(); for my $i (0..$OFF+1) { my $b = R(1); for my $q (0..$i-1) { $b = $b * ($tau + $q) / ($q + 1); } $z0 = ladd($z0, lmon(1-$i, $b)); }
# Newton iteration for z = z_0 + delta (bivariate): F(z) := alpha Y_0(z/alpha) - sum_j P_j z^j - lambda - tau = 0
# Y_0(z) = z - sigma sum_k Pi_k z^{-k}/(k+1).  Need z^{-k}: z = lambda (1 + ...) -> invert
my $alB = bz(); for my $m (0..$M) { $alB->[$m] = lmon(0, $W[$m]); }
my $z = bz(); $z->[0] = [ @$z0 ];
sub binv_lead_lambda { my $zz = shift; # 1/zz for zz = lambda*(1 + small)
  my $q = bz(); for my $m (0..$M) { $q->[$m] = lmul($zz->[$m], lmon(-1)); } my $d = badd($q, bconst(R(1)), R(-1));
  my $inv = bconst(R(1)); my $pw = bconst(R(1)); for my $k (1..$OFF + $TOPP + 2) { $pw = bmul($pw, $d); $inv = badd($inv, $pw, R((-1)**$k)); }
  my $r = bz(); for my $m (0..$M) { $r->[$m] = lmul($inv->[$m], lmon(-1)); } return $r; }
for my $it (1..$M+2) {
  my $zi = binv_lead_lambda($z);
  my $Fz = [ map { [ @$_ ] } @$z ];   # alpha Y_0(z/alpha) = z - sigma sum_k Pi_k alpha^{k+1} z^{-k}/(k+1)
  my $zk = bconst(R(1)); my $alp = [ map { [ @$_ ] } @$alB ];
  for my $k (1..$OFF + 2) { $zk = bmul($zk, $zi); $alp = bmul($alp, $alB); $Fz = badd($Fz, bmul($alp, $zk), -$sig*Pi_at($k,$tau)/($k+1)); }
  my $zj = [ map { [ @$_ ] } @$z ]; for my $j (2..$M+1) { $zj = bmul($zj, $z); my $Pb = bz(); for my $m (0..$M) { $Pb->[$m] = lmon(0, $Pj{$j}[$m]); } $Fz = badd($Fz, bmul($Pb, $zj), R(-1)); }
  $Fz = badd($Fz, bconst($tau), R(-1)); $Fz->[0] = ladd($Fz->[0], lmon(1, R(-1)));
  # derivative ~ 1 at leading order: Newton with approximate derivative Y_0'(z_0) ~ 1 + O(lambda^-2); use simple fixed-point z <- z - F(z)/Y0'(z0)
  my $dY = lmon(0, R(1)); # Y_0'(z) = 1 + sigma sum_k k Pi_k z^{-k-1}/(k+1): evaluate on z_0
  my $zi0 = binv_lead_lambda([ $z0, map { lz() } 1..$M ])->[0]; my $p = lmon(0, R(1));
  for my $k (1..$OFF + 2) { $p = lmul($p, $zi0); my $t2 = lmul($p, $zi0); $dY = ladd($dY, $t2, $sig*Pi_at($k,$tau)*$k/($k+1)); }
  # 1/dY
  my $dd = ladd($dY, lmon(0, R(1)), R(-1)); my $idY = lmon(0, R(1)); my $pp = lmon(0, R(1)); for my $k (1..$OFF+$TOPP) { $pp = lmul($pp, $dd); $idY = ladd($idY, $pp, R((-1)**$k)); }
  my $corr = bz(); for my $m (0..$M) { $corr->[$m] = lmul($Fz->[$m], $idY); } $z = badd($z, $corr, R(-1)); }
# K = log(z u) = log(z/z_0): z/z_0 = 1 + small
my $zi0b = binv_lead_lambda([ $z0, map { lz() } 1..$M ]);
my $ratio = bz(); for my $m (0..$M) { $ratio->[$m] = lmul($z->[$m], $zi0b->[0]); }
my $sm = badd($ratio, bconst(R(1)), R(-1)); my $K = bz(); my $pw = bconst(R(1));
for my $k (1..$M) { $pw = bmul($pw, $sm); $K = badd($K, $pw, R((-1)**($k+1))/$k); }
my @Bn = (R(1)); for my $m (1..40) { my $s = R(0); for my $k (0..$m-1) { $s = $s + binom($m+1, $k) * $Bn[$k]; } $Bn[$m] = -$s/($m+1); }
my $bad = 0;
for my $k (1..$M) { my $kap = $k == 1 ? $tau + R('1/2') : $Bn[$k]/($k*fact($k));
  my @pred = map { R(0) } 0..$k; for my $i (0..$k) { $pred[$i] = $kap * ( ($i == $k ? 1+$tau : R(0)) - $tau * binom($k, $i) * R((-1)**($k-$i)) ); }
  my @got = map { $K->[$k][$_ + $OFF] } 0..$k;
  my $ok = 1; for my $i (0..$k) { $ok = 0 if $got[$i] != $pred[$i]; } $bad++ unless $ok;
  print "t^$k polynomial part: got (", join(', ', map { "$_" } @got), ")  T1 (", join(', ', map { "$_" } @pred), ")  ", ($ok ? "OK" : "DIFF"), "\n";
  my @hi = map { $K->[$k][$_ + $OFF] } ($k+1 .. $TOPP); print "   higher powers present: ", (grep { !$_->is_zero } @hi) ? "yes" : "no", "\n"; }
print "t_ansatz: ", ($bad ? "$bad orders differ" : "all orders agree"), "\n";
