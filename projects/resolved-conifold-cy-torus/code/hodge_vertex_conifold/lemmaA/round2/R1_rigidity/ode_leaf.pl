#!/usr/bin/perl
# ode_leaf.pl -- independent, polynomial-time computation of the leaf-curve functional Phi(t) = F_G(P(t)) of a genus-G
# linearised KP solution (Route A, Theorem D, at every point of the slice P_1 = 0), by integrating along the leaf curve:
#   U_n(t) = (d_{P_n} F)(P(t)) (n >= 2),  U_1 = (E F)(P(t)),  G_{mn} = two-point function at P(t),
#   Hirota at P(t):  sum G_{mn} z^-m w^-n = (V(z) - V(w))/(Y_t(w) - Y_t(z)),  Y_t(z) = alpha(t) Y_0(z/alpha(t))   [genus 0 on the slice],
#   string on the slice:  V_n = G_{1n} = n U_n + sum_{j>=2} j P_j G_{jn},
#   d/dt U_n = sum_{i>=2} P_i'(t) G_{in},   d/dt Phi = sum_{i>=2} P_i'(t) U_i,   U_n(0) = Y(n) = y_{n-1}.
# In the rescaled y-coordinate (u^ = alpha/z = y(1-y)^tau) the two-point function is Hankel:
#   G_{ab} = alpha^{a+b-1} sum_{i<=b, j<=a} c_{a,j} c_{b,i} F_{i+j-1},  F_k = sum_n alpha^{-n} d_{k,n} V_n,
#   c_{n,i} = [u^n] y^i = (i/n) binom(n tau + n - i - 1, n - i),  d_{k,n} = [y^k] (y(1-y)^tau)^n = (-1)^{k-n} binom(n tau, k-n).
# Compares Phi(t) with the closed form (T1):  [t^k] Phi = kappa_k sum_L y_L (-1)^{k-L} [ (1+tau) binom((L+1)tau, k-L) - tau binom((L+1)(1+tau)+k-L-1, k-L) ].
# Modular arithmetic (p = 2^31-1).  Usage: perl ode_leaf.pl K tau(a/b) [seed]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod;
my ($K, $ts, $seed) = @ARGV; $seed //= 3; my ($ta, $tb) = split m{/}, $ts; $tb //= 1; my $tau = mfrac($ta, $tb); my $sig = mm($tau, md($tau+1));
my $NMAX = $K + 2;
my @fct = (1); $fct[$_] = mm($fct[$_-1], $_) for 1..4*$K+10;
sub ifct { mi($fct[$_[0]]) }
sub gbin { my ($x, $m) = @_; return 0 if $m < 0; my $v = 1; for my $i (0..$m-1) { $v = mm($v, md($x - $i)); } return mm($v, ifct($m)); }
sub Pi { my $j = shift; return 1 if $j <= 1; my $v = 1; for my $b (1..$j-1) { $v = mm($v, md($j*$tau + $b)); } return mm($v, ifct($j-1)); }
# series helpers (arrays mod p, length K+1)
sub smul { my ($a, $b) = @_; my @c = (0) x ($K+1); for my $i (0..$K) { next unless $a->[$i]; for my $j (0..$K-$i) { $c[$i+$j] = ($c[$i+$j] + mm($a->[$i], $b->[$j])) % $P if $b->[$j]; } } return \@c; }
sub sinv { my $a = shift; my @b = (mi($a->[0])); for my $n (1..$K) { my $s = 0; for my $k (1..$n) { $s = ($s + mm($a->[$k], $b[$n-$k])) % $P; } $b[$n] = md(-mm($s, $b[0])); } return \@b; }
sub scomp { my ($f, $g) = @_; my @res = (0) x ($K+1); my $pw = [1, (0) x $K]; for my $k (0..$K) { for my $i (0..$K) { $res[$i] = ($res[$i] + mm($f->[$k] // 0, $pw->[$i])) % $P; } $pw = smul($pw, $g); } return \@res; }
sub hyp { my ($a, $b, $c) = @_; my @s = (1); my $t = 1; for my $n (1..$K) { $t = mm($t, mm(mm(md($a + $n - 1), md($b + $n - 1)), mi(mm(md($c + $n - 1), $n)))); push @s, $t; } return \@s; }
# leaf curve
my @Tt = (0); for my $n (1..$K) { push @Tt, mm(md(($n % 2) ? 1 : -1), ifct($n)); }
my $h1 = scomp(hyp(md(-$tau), md(-$tau), 1), \@Tt);
my %Pj; my %dPj;
for my $j (2..$NMAX) { my $hj = scomp(hyp(md(-$j*$tau), md(-$j*$tau), $j), \@Tt); my $p = [1, (0) x $K]; $p = smul($p, $h1) for 1..$j; my $rho = smul($hj, sinv($p));
  my $cc = mm(md(($j % 2) ? -1 : 1), mm(Pi($j), mi($j))); my @s = (0) x ($K+1); for my $m ($j-1..$K) { $s[$m] = mm($cc, $rho->[$m-$j+1]); } $Pj{$j} = \@s;
  my @d = (0) x ($K+1); for my $m (0..$K-1) { $d[$m] = mm($m+1, $s[$m+1]); } $dPj{$j} = \@d; }
# alpha(t): W = 1 + sum_j Pi_j P_j W^j
my @al = (1, (0) x $K);
for my $it (0..$K+1) { my @nw = (1, (0) x $K); my $wp = [@al]; for my $j (2..$NMAX) { $wp = smul($wp, \@al); my $q = smul([ map { mm(Pi($j), $_) } @{$Pj{$j}} ], $wp); $nw[$_] = ($nw[$_] + $q->[$_]) % $P for 0..$K; } @al = @nw; }
my @alp = ([1, (0) x $K]); for my $e (1..2*$NMAX+2) { $alp[$e] = smul($alp[$e-1], \@al); }   # alpha^e
my $ali = sinv(\@al); my @alm = ([1, (0) x $K]); for my $e (1..2*$NMAX+2) { $alm[$e] = smul($alm[$e-1], $ali); }   # alpha^{-e}
# Lagrange coefficients
my $KM = 2*$NMAX + 2;
my %c; for my $n (1..$KM) { for my $i (1..$n) { $c{"$n,$i"} = mm(mm($i, mi($n)), gbin(md($n*$tau + $n - $i - 1), $n - $i)); } }
my %d; for my $k (1..$KM) { for my $n (1..$k) { $d{"$k,$n"} = mm(md((($k-$n) % 2) ? -1 : 1), gbin(mm($n, $tau), $k - $n)); } }
# one-leaf data
srand($seed); my @y = (0); push @y, int(rand(2000)) + 1 for 1..$KM;   # y_L, L = 1..KM ; U_n(0) = y_{n-1}
# state: U[n][m], V[n][m], F[k][m], G{a,b}[m]  with level truncation (n-1)+m <= K etc.
my (@U, @V, @F, %G); my @Phi = (0) x ($K+2);
for my $n (1..$NMAX) { $U[$n] = [ (0) x ($K+1) ]; $V[$n] = [ (0) x ($K+1) ]; } for my $k (1..$KM) { $F[$k] = [ (0) x ($K+1) ]; }
for my $n (2..$NMAX) { $U[$n][0] = $y[$n-1]; } $U[1][0] = 0;
for my $m (0..$K) {
  # V_n^{(m)} = n U_n^{(m)} + sum_{j>=2} j sum_{m'>=1} P_j^{(m')} G_{jn}^{(m-m')}
  for my $n (1..$NMAX) { next if ($n - 1) + $m > $K + 1; my $v = mm($n, $U[$n][$m]);
    for my $j (2..$NMAX) { for my $mp ($j-1..$m) { my $pj = $Pj{$j}[$mp]; next unless $pj; my $g = $G{"$j,$n"}; next unless $g; $v = ($v + mm(mm($j, $pj), $g->[$m-$mp])) % $P; } }
    $V[$n][$m] = $v; }
  # F_k^{(m)} = sum_n sum_{m'} (alpha^{-n})^{(m')} d_{k,n} V_n^{(m-m')}
  for my $k (1..$KM) { next if ($k - 1) + $m > 2*$K + 2; my $s = 0;
    for my $n (1..($k < $NMAX ? $k : $NMAX)) { my $dk = $d{"$k,$n"}; next unless $dk; for my $mp (0..$m) { my $a = $alm[$n][$mp]; next unless $a; $s = ($s + mm(mm($a, $dk), $V[$n][$m-$mp])) % $P if $V[$n][$m-$mp]; } }
    $F[$k][$m] = $s; }
  # G_{ab}^{(m)} = sum_{m'} (alpha^{a+b-1})^{(m')} Ghat_{ab}^{(m-m')},  Ghat_{ab} = sum_{i<=b, j<=a} c_{a,j} c_{b,i} F_{i+j-1}
  for my $a (1..$NMAX) { for my $b (1..$NMAX) { next if ($a-1)+($b-1)+$m > $K + 1;
      my $gh = 0; for my $j (1..$a) { my $ca = $c{"$a,$j"}; next unless $ca; for my $i (1..$b) { my $cb = $c{"$b,$i"}; next unless $cb; $gh = ($gh + mm(mm($ca, $cb), $F[$i+$j-1][$m])) % $P; } }
      $G{"$a,$b"}[$m] = 0 unless defined $G{"$a,$b"}[$m]; $G{"$a,$b"}[$m] = $gh; } }   # store Ghat first, multiply by alpha powers below
  # convert Ghat -> G at order m: G^{(m)} = sum_{m'} alpha^{a+b-1}(m') Ghat^{(m-m')}
  for my $a (1..$NMAX) { for my $b (1..$NMAX) { next if ($a-1)+($b-1)+$m > $K + 1; my $gh = $G{"$a,$b"}; my $hat = $G{"H$a,$b"} //= [];
      $hat->[$m] = $gh->[$m]; my $s = 0; my $e = $a + $b - 1; for my $mp (0..$m) { my $ap = $alp[$e][$mp]; next unless $ap; $s = ($s + mm($ap, $hat->[$m-$mp] // 0)) % $P; } $gh->[$m] = $s; } }
  last if $m == $K;
  # U_n^{(m+1)} = (1/(m+1)) sum_i sum_{m'} P_i'^{(m')} G_{in}^{(m-m')}
  for my $n (1..$NMAX) { next if ($n - 1) + $m + 1 > $K + 1; my $s = 0;
    for my $i (2..$NMAX) { for my $mp ($i-2..$m) { my $dp = $dPj{$i}[$mp]; next unless $dp; my $g = $G{"$i,$n"}; next unless $g && defined $g->[$m-$mp]; $s = ($s + mm($dp, $g->[$m-$mp])) % $P; } }
    $U[$n][$m+1] = mm($s, mi($m+1)); }
  # Phi^{(m+1)} = (1/(m+1)) sum_i sum_{m'} P_i'^{(m')} U_i^{(m-m')}
  my $s = 0; for my $i (2..$NMAX) { for my $mp ($i-2..$m) { my $dp = $dPj{$i}[$mp]; next unless $dp; $s = ($s + mm($dp, $U[$i][$m-$mp])) % $P; } }
  $Phi[$m+1] = mm($s, mi($m+1));
}
# closed form T1
my @B = (1); for my $mm (1..$K) { my $s = 0; for my $k (0..$mm-1) { $s = ($s + mm(mm($fct[$mm+1], mm(ifct($k), ifct($mm+1-$k))), $B[$k])) % $P; } $B[$mm] = md(-mm($s, mi($mm+1))); }
my $bad = 0;
for my $k (1..$K) { my $kap = $k == 1 ? md($tau + mi(2)) : mm($B[$k], mi(mm($k, $fct[$k])));
  my $pred = 0; for my $L (1..$k) { my $N = $L + 1; my $dd = $k - $L;
    my $g1 = gbin(mm($N, $tau), $dd); my $g2 = gbin(md(mm($N, md(1+$tau)) + $dd - 1), $dd);
    my $bt = mm(md(($dd % 2) ? -1 : 1), md(mm(md(1+$tau), $g1) - mm($tau, $g2))); $pred = ($pred + mm($bt, $y[$L])) % $P; }
  $pred = mm($pred, $kap); $bad++ if $pred != $Phi[$k]; printf "t^%2d: ODE %10d  closed form %10d  %s\n", $k, $Phi[$k], $pred, ($pred == $Phi[$k] ? "OK" : "DIFF") if $k <= 6 || $pred != $Phi[$k]; }
print "ode_leaf K=$K tau=$ts: ", ($bad ? "$bad mismatches" : "closed form T1 confirmed to order $K"), "\n";
