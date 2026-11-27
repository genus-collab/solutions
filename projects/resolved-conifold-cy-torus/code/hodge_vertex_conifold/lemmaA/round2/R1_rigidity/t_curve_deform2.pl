#!/usr/bin/perl
# t_curve_deform2.pl -- is lambda + tau = Y_t(z) - sum_{j>=1} c_j(t) z^j on the curve z = z_t(lambda) := 1/u_t(lambda)?
#   u_t(lambda) = t (e^{t(lambda-1)} - 1)^tau (e^{t lambda} - 1)^{-1-tau} e^{-tau t (lambda + tau)},
#   Y_t(z) = alpha(t) Y_0(z/alpha(t)),  Y_0(z) = z - sigma sum_k Pi_k z^{-k}/(k+1),  alpha = genus-0 scalar on the leaf curve.
# Double series: t-order m = 0..M; at each order a Laurent polynomial/series in lambda, truncated at lambda^{-NEG}.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M, $NEG) = @ARGV; $taus //= '2/7'; $M //= 3; $NEG //= 8; my $tau = R($taus); my $sig = $tau*(1+$tau);
my $TOP = $M + 2;   # max positive power of lambda tracked
# Laurent series in lambda: hash power => coeff, powers from -NEG to TOP*? ; use arrays offset
my $OFF = $NEG; my $LEN = $NEG + 2*$M + 6;
sub lz { [ map { R(0) } 0..$LEN ] }
sub lmul { my ($a, $b) = @_; my $c = lz(); for my $i (0..$LEN) { next if $a->[$i]->is_zero; for my $j (0..$LEN) { next if $b->[$j]->is_zero; my $k = $i + $j - $OFF; next if $k < 0 || $k > $LEN; $c->[$k] = $c->[$k] + $a->[$i]*$b->[$j]; } } return $c; }
sub ladd { my ($a, $b, $s) = @_; $s //= R(1); return [ map { $a->[$_] + $s*$b->[$_] } 0..$LEN ]; }
sub lmon { my ($p, $co) = @_; my $v = lz(); $v->[$p + $OFF] = $co // R(1) if $p + $OFF >= 0 && $p + $OFF <= $LEN; return $v; }
# bivariate: array over t-order of Laurent series
sub bz { [ map { lz() } 0..$M ] }
sub bmul { my ($a, $b) = @_; my $c = bz(); for my $i (0..$M) { for my $j (0..$M-$i) { $c->[$i+$j] = ladd($c->[$i+$j], lmul($a->[$i], $b->[$j])); } } return $c; }
sub badd { my ($a, $b, $s) = @_; return [ map { ladd($a->[$_], $b->[$_], $s) } 0..$M ]; }
sub bconst { my $v = bz(); $v->[0] = lmon(0, $_[0]); return $v; }
# power of a bivariate series of the form X0 * (1 + small) with X0 a Laurent series having leading term: use log/exp on the "1+small" factor
# z_t(lambda) = (e^{t lambda}-1)^{1+tau} (e^{t(lambda-1)}-1)^{-tau} e^{tau t (lambda+tau)} / t
#   = [lambda (1 + E1)]^{1+tau} [(lambda-1)(1 + E2)]^{-tau} e^{tau t (lambda + tau)},  E1 = (e^{t lambda}-1)/(t lambda) - 1, E2 likewise with lambda-1
# log z_t = log z_0(lambda) + (1+tau) log(1+E1) - tau log(1+E2) + tau t (lambda + tau),   z_0 = lambda^{1+tau}(lambda-1)^{-tau} = lambda (1 - 1/lambda)^{-tau}
# (1+tau)log((e^{x}-1)/x) with x = t lambda : series sum_n K_n x^n, K_n = [x^n] log((e^x-1)/x)
my @Bn = (R(1)); for my $m (1..40) { my $s = R(0); for my $k (0..$m-1) { $s = $s + binom($m+1, $k) * $Bn[$k]; } $Bn[$m] = -$s/($m+1); }
my @Kn = (R(0)); for my $n (1..$M+2) { $Kn[$n] = $n == 1 ? R('1/2') : $Bn[$n]/($n*fact($n)); }
my $logratio = bz();   # log(z_t/z_0)
for my $n (1..$M) { # t^n coefficient: (1+tau) K_n lambda^n - tau K_n (lambda-1)^n + [n==1] tau (lambda + tau)
  my $v = ladd(lmon($n, (1+$tau)*$Kn[$n]), lz());
  for my $i (0..$n) { $v = ladd($v, lmon($i, -$tau*$Kn[$n]*binom($n,$i)*R((-1)**($n-$i)))); }
  if ($n == 1) { $v = ladd($v, lmon(1, $tau)); $v = ladd($v, lmon(0, $tau*$tau)); }
  $logratio->[$n] = $v; }
# exp of logratio (no t^0 term)
my $ex = bconst(R(1)); my $term = bconst(R(1));
for my $k (1..$M) { $term = bmul($term, $logratio); my $f = R(1)/fact($k); $ex = badd($ex, $term, $f); }
# z_0(lambda) = lambda (1 - 1/lambda)^{-tau} = sum_i binom(tau+i-1, i) lambda^{1-i}
my $z0 = lz(); for my $i (0..$OFF+1) { my $b = R(1); for my $q (0..$i-1) { $b = $b * ($tau + $q) / ($q + 1); } $z0 = ladd($z0, lmon(1-$i, $b)); }
my $zt = bz(); for my $m (0..$M) { $zt->[$m] = lmul($ex->[$m], $z0); }
# alpha(t): genus-0 scalar (from t_curve_deform.pl logic, recomputed here)
my $x = -$tau/(1+$tau); my $c = 1 + $x; my $N = $M + 2; my $Tt = T_of_t_series($N+1);
my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $N+1), $Tt, $N+1);
my @u = (R(0), 4*$x/($c*$c), map { R(0) } 2..$N+1); my @omu = map { -$_ } @u; $omu[0] = R(1);
my $sq = series_pow(\@omu, R('1/2'), $N+1); my @den = map { $c * $_ } @$sq; $den[0] = $den[0] + $c;
my $ss = series_inv(\@den, $N+1); $ss = [ map { 2*$_ } @$ss ];
my @Tsx = (R(0), @$ss[0..$N]); my @omTs = map { -$_ } @Tsx; $omTs[0] = R(1);
my $dd = series_mul($ss, series_inv(\@omTs, $N+1), $N+1); my @dphi = (R(0)); $dphi[$_] = $dd->[$_-1]/$_ for 1..$N+1;
my $e = series_exp([ map { $_ / $c } @dphi ], $N+1); my @zT = (R(0), @$e[0..$N]); my $zeta = series_compose(\@zT, $Tt, $N+1);
my $al = series_mul([ map { $zeta->[$_+1] } 0..$N ], $h1, $N);
# Y_t(z_t) = z_t - sigma sum_k Pi_k alpha^{k+1} z_t^{-k}/(k+1)
# z_t^{-1}: z_t = lambda * (1 + ...) ; invert bivariate series
my $lead = lz(); $lead = lmon(1, R(1));   # lambda
my $zt_over = bz(); for my $m (0..$M) { $zt_over->[$m] = lmul($zt->[$m], lmon(-1)); }   # z_t/lambda = 1 + O(1/lambda) + O(t)
my $one = bconst(R(1)); my $dlt = badd($zt_over, $one, R(-1));
my $inv = bconst(R(1)); my $pw = bconst(R(1)); for my $k (1..$NEG + 2*$M + 2) { $pw = bmul($pw, $dlt); $inv = badd($inv, $pw, R((-1)**$k)); }
my $zinv = bz(); for my $m (0..$M) { $zinv->[$m] = lmul($inv->[$m], lmon(-1)); }   # 1/z_t
my $Yt = [ map { [ @$_ ] } @$zt ];
my $zk = bconst(R(1)); my $alser = bz(); for my $m (0..$M) { $alser->[$m] = lmon(0, $al->[$m]); }
my $alpow = [ map { [ @$_ ] } @$alser ];   # alpha^{k+1}, start with alpha^1 then multiply
for my $k (1..$NEG + 2) { $zk = bmul($zk, $zinv); $alpow = bmul($alpow, $alser); my $f = -$sig * Pi_at($k, $tau) / ($k + 1); $Yt = badd($Yt, bmul($alpow, $zk), $f); }
# Delta = Y_t(z_t) - tau - lambda
my $Delta = badd($Yt, bconst($tau), R(-1)); $Delta->[0] = ladd($Delta->[0], lmon(1, R(-1)));
# fit Delta = sum_j c_j(t) z_t^j : at each t-order, peel the highest powers of lambda
my @ztp = (bconst(R(1))); for my $j (1..$M+1) { $ztp[$j] = bmul($ztp[$j-1], $zt); }
my %cj; my $res = [ map { [ @$_ ] } @$Delta ];
for my $m (0..$M) { for my $p (reverse 0..$M+1) {   # power of lambda p at t-order m
    my $co = $res->[$m][$p + $OFF]; next if $co->is_zero; # remove with c_{p,m} z_t^p (leading lambda^p at t^0)
    $cj{"$p,$m"} = $co; my $sh = bz(); for my $mm (0..$M-$m) { $sh->[$m+$mm] = $ztp[$p][$mm]; } $res = badd($res, $sh, -$co); } }
for my $j (0..$M+1) { my @c = map { $cj{"$j,$_"} // R(0) } 0..$M; print "c_$j(t) = ", join(', ', map { "$_" } @c), "\n"; }
my $maxres = 0; for my $m (0..$M) { for my $p (-$NEG+$M+3 .. $M+2) { my $i = $p + $OFF; next if $i < 0 || $i > $LEN; $maxres++ unless $res->[$m][$i]->is_zero; } }
print "nonzero residual coefficients (powers above the truncation zone): $maxres\n";
for my $m (0..$M) { print "  residual t^$m: ", join(' ', map { my $i = $_ + $OFF; ($i >= 0 && $i <= $LEN && !$res->[$m][$i]->is_zero) ? "l^$_:" . $res->[$m][$i] : () } (-4..2)), "\n"; }
print "Delta(lambda,t) = Y_t(z_t(lambda)) - lambda - tau, coefficients (powers of lambda -4..", $M+2, "):\n";
for my $m (0..$M) { print "  t^$m: ", join(' ', map { my $i = $_ + $OFF; ($i >= 0 && $i <= $LEN && !$Delta->[$m][$i]->is_zero) ? "l^$_:" . $Delta->[$m][$i] : () } (-4..$M+2)), "\n"; }
