#!/usr/bin/perl
# xy_onepoint_generic.pl tau G jmax : n = 1 x-y swap formula for ARBITRARY rational tau (no finite products):
#   [x^j] F_{g,1} = (1/j^2) [hb^{2g}] [v^{-1}]_{v=oo} exp(-j S(j hb d)^{-1}... )  realised as
#   Q_j(v) = x^{-j} exp(-j (Sc(j hb d) xhat - X)),  xhat = tau S(hb d/tau) log v - e2 S(hb d/e2) log(1+v),  Sc(t) = 1/S(t) = 2 sinh(t/2)/t,
# expanded at v = oo: log terms enter only through d^{2k} log(v - a) = -(2k-1)! (v - a)^{-2k}.  x^{-j} = v^{j} (1 + 1/v)^{e2 j}.
# Compared with (6.1) through the dictionary P2: sigma^g Xi_g(j) = sigma^{2g-1} [x^j]F_{g,1} / Pi_j.
use strict; use warnings; use Math::BigRat;
sub Q { Math::BigRat->new(@_) }
my ($tau, $G, $jmax) = @ARGV; $tau = Q($tau // "1/2"); $G //= 3; $jmax //= 6;
my $e2 = 1 + $tau; my $sig = $tau * $e2; my $e1 = -$tau;
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; return $r; }
sub bern { my $n = shift; my @B = (Q(1)); for my $m (1..$n) { my $s = Q(0); for my $k (0..$m-1) { $s += fact($m+1)/(fact($k)*fact($m+1-$k)) * $B[$k]; } $B[$m] = -$s/($m+1); } return $B[$n]; }
my $K = 2*$G;           # hb order
# S(t) = sum s_k t^{2k}, Sc(t) = 1/S(t) = sum c_k t^{2k}
my @s = map { $_ % 2 ? Q(0) : Q(0) } 0..$K;
for my $k (0..$G) { $s[2*$k] = bern(2*$k) * (Q(1)/Q(2))**0; }
# S(t) = (t/2)/sinh(t/2) = sum_k B_{2k}(1/2) t^{2k}/(2k)!,  B_{2k}(1/2) = (2^{1-2k} - 1) B_{2k}
for my $k (0..$G) { $s[2*$k] = (Q(2)**(1-2*$k) - 1) * bern(2*$k) / fact(2*$k); }
my @sc = (Q(1)); for my $n (1..$K) { my $acc = Q(0); for my $m (1..$n) { $acc += $s[$m] * $sc[$n-$m]; } $sc[$n] = -$acc; }
# Laurent series at v = oo in w = 1/v, with coefficients polynomial in hb: hash  "pw,ph" => coef  (pw = power of w, may be negative)
sub smul { my ($A, $B, $wmax) = @_; my %C;
  for my $ka (keys %$A) { my ($a1, $a2) = split /,/, $ka; for my $kb (keys %$B) { my ($b1, $b2) = split /,/, $kb;
      my ($c1, $c2) = ($a1 + $b1, $a2 + $b2); next if $c2 > $K || $c1 > $wmax; my $k = "$c1,$c2"; $C{$k} = ($C{$k} // Q(0)) + $A->{$ka} * $B->{$kb}; } }
  for (keys %C) { delete $C{$_} if $C{$_}->is_zero; } return \%C; }
sub sexp { my ($A, $wmax) = @_; my %R = ("0,0" => Q(1)); my %T = ("0,0" => Q(1));
  for my $n (1..(2*$K + $wmax + 4)) { my $Tn = smul(\%T, $A, $wmax); %T = map { $_ => $Tn->{$_} / $n } keys %$Tn; last unless %T;
    for (keys %T) { $R{$_} = ($R{$_} // Q(0)) + $T{$_}; } }
  return \%R; }
# d^{2k} log(v - a) at v = oo:  -(2k-1)! (v - a)^{-2k} = -(2k-1)! w^{2k} (1 - a w)^{-2k}
sub dlog { my ($kk, $a, $wmax) = @_; my %R; my $c = -fact(2*$kk - 1);
  for my $m (0..($wmax - 2*$kk)) { my $bin = fact(2*$kk + $m - 1) / (fact($m) * fact(2*$kk - 1)); $R{(2*$kk + $m) . ",0"} = $c * $bin * $a**$m; } return \%R; }
for my $j (1..$jmax) {
  my $wmax = $j + 2;
  # exponent E = -j (Sc(j hb d) xhat - X): only terms with total d-order 2k >= 2 survive
  # Sc(j hb d) xhat = sum_{a,b} sc_a (j hb)^a d^a [ tau s_b (hb/tau)^b d^b log v - e2 s_b (hb/e2)^b d^b log(1+v) ]
  my %E;
  for my $a (0..$K) { next if $sc[$a]->is_zero; for my $b (0..$K - $a) { next if $s[$b]->is_zero; next if $a + $b == 0; next if ($a + $b) % 2;
      my $kk = ($a + $b) / 2;
      my $coef1 = $sc[$a] * Q($j)**$a * $tau * $s[$b] * $tau**(-$b);
      my $coef2 = -$sc[$a] * Q($j)**$a * $e2 * $s[$b] * $e2**(-$b);
      my $L1 = dlog($kk, Q(0), $wmax); my $L2 = dlog($kk, Q(-1), $wmax);
      for my $key (keys %$L1) { my ($pw) = split /,/, $key; my $k2 = "$pw," . ($a + $b); $E{$k2} = ($E{$k2} // Q(0)) - $j * $coef1 * $L1->{$key}; }
      for my $key (keys %$L2) { my ($pw) = split /,/, $key; my $k2 = "$pw," . ($a + $b); $E{$k2} = ($E{$k2} // Q(0)) - $j * $coef2 * $L2->{$key}; } } }
  for (keys %E) { delete $E{$_} if $E{$_}->is_zero; }
  my $EX = sexp(\%E, $wmax);
  # x^{-j} = v^j (1 + w)^{e2 j}:  [v^{-1}] Q = [w^{1}] of w^{-j} (1+w)^{e2 j} EX  = [w^{j+1}] (1+w)^{e2 j} EX
  my %P; my $bin = Q(1); for my $m (0..$wmax) { $P{"$m,0"} = $bin; $bin = $bin * ($e2 * $j - $m) / ($m + 1); }
  my $Prod = smul(\%P, $EX, $wmax);
  my $Pi = Q(1); $Pi *= ($j * $tau + $_) for 1..$j-1; $Pi /= fact($j-1);
  for my $g (1..$G) { my $val = ($Prod->{($j+1) . "," . (2*$g)} // Q(0)) / ($j * $j);
    my $lhs = $sig**(2*$g - 1) * $val / $Pi;
    # (6.1)
    my @lc = map { Q(0) } 0..$K; for my $n (1..$G) { my $cn = -bern(2*$n) / (2*$n * fact(2*$n)); $lc[2*$n] = $cn * ((1 + $j*$e2) * $e1**(2*$n) + (1 + $j*$e1) * $e2**(2*$n)); }
    my @ex = (Q(1)); for my $n (1..$K) { my $acc = Q(0); for my $m (1..$n) { $acc += $m * $lc[$m] * $ex[$n-$m]; } $ex[$n] = $acc / $n; }
    my $pre = Q(1); $pre *= ($j - $_) for 1..(2*$g - 2); my $rhs = $pre * $ex[2*$g];
    printf "tau=%s j=%d g=%d: xy %s  (6.1) %s  %s\n", $tau, $j, $g, $lhs, $rhs, ($lhs == $rhs ? "OK" : "MISMATCH"); } }
