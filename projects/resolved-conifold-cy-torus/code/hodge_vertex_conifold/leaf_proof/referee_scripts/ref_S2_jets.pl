#!/usr/bin/perl
# Referee check for S2 (independent of S1/S2 closed forms): build the code's rows exactly as rc_lib.m/hvc_lib.m do
# (Wick series in y = log U at a critical point, frame-B twist h_w^-1, code's literal Delta = 1/(e_1^T eta e_1),
# H-row recursion), as exact JETS in e = t - t0 at a generic point t0 != 0 where everything is rational
# (parametrise by s0 = s*(t0) rational).  Then test:
#   (V1) Delta_i = x^2 W2(U_i);
#   (V5) second QDE row at BOTH idempotents: (D_i r_{1,i})[m-1] = L_j (r_{1,i}[m] - L_i r_{0,i}[m]);
#   (V6) symplecticity / orthogonality: Sum_{p+q=n} (-1)^q <v_i[p], v_j[q]> = delta_{n0} delta_ij Delta_i.
# Usage: perl ref_S2_jets.pl [M] [J]   (M = max z-order, J = jet order)
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $M = @ARGV ? $ARGV[0] : 6;
my $J = @ARGV > 1 ? $ARGV[1] : 4;
my ($ntest, $nfail) = (0, 0);
sub ok { my ($c, $msg) = @_; $ntest++; unless ($c) { $nfail++; print "FAIL: $msg\n"; } }

# ---------------- jets: arrays [a0..aJ], a_k = coefficient of e^k
sub jz { [ map { R(0) } 0..$J ] }
sub jc { my $a = jz(); $a->[0] = R($_[0]); $a }
sub jadd { my ($a,$b)=@_; [ map { $a->[$_] + $b->[$_] } 0..$J ] }
sub jsub { my ($a,$b)=@_; [ map { $a->[$_] - $b->[$_] } 0..$J ] }
sub jsc  { my ($s,$a)=@_; [ map { $a->[$_] * $s } 0..$J ] }
sub jmul { my ($a,$b)=@_; my $c = jz(); for my $i (0..$J) { next if $a->[$i]->is_zero; for my $k (0..$J-$i) { $c->[$i+$k] += $a->[$i]*$b->[$k]; } } $c }
sub jinv { my $a = shift; die "jinv 0" if $a->[0]->is_zero; my $b = jz(); $b->[0] = 1/$a->[0];
  for my $n (1..$J) { my $s = R(0); $s += $a->[$_]*$b->[$n-$_] for 1..$n; $b->[$n] = -$s/$a->[0]; } $b }
sub jder { my $a = shift; my $b = jz(); $b->[$_-1] = $a->[$_]*$_ for 1..$J; $b }   # top coefficient becomes invalid
sub jexp { my $a = shift; my $e = jz(); $e->[0] = R(1); my $a0 = $a->[0]; die "jexp needs a0=0" unless $a0->is_zero;
  for my $n (1..$J) { my $s = R(0); $s += $a->[$_]*$e->[$n-$_]*$_ for 1..$n; $e->[$n] = $s/$n; } $e }
sub jeq { my ($a,$b,$upto) = @_; for my $i (0..$upto) { return 0 unless $a->[$i] == $b->[$i]; } 1 }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }

# integer partitions of N into exactly k parts, each in [lo, hi], as lists (non-increasing)
sub parts { my ($N, $k, $lo, $hi) = @_; my @out;
  return ([]) if $k == 0 && $N == 0; return () if $k == 0;
  for my $p (reverse $lo..($hi < $N ? $hi : $N)) { next if $p*$k < $N; next if $lo*($k-1) > $N - $p;
    for my $rest (parts($N - $p, $k - 1, $lo, $p)) { push @out, [$p, @$rest]; } }
  @out }

my @bern = (undef, R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"));   # B_2k

for my $pt ([R(2), R("2/5")], [R("3/5"), R("3/4")], [R(-3), R("-2/5")], [R("1/3"), R("6/7")]) {
  my ($x, $s0) = @$pt; my $c = 1 + $x;
  my $tau0 = ($c*$s0 - 1)/($x*$s0*$s0);
  # tau(t0 + e) = 1 - (1 - tau0) e^{-e};  Q = Q0 e^{e}
  my $em = jz(); $em->[$_] = R((($_ % 2) ? -1 : 1))/fact($_) for 0..$J;     # e^{-e}
  my $ep = jz(); $ep->[$_] = R(1)/fact($_) for 0..$J;                       # e^{e}
  my $tau = jsub(jc(1), jsc(1 - $tau0, $em));
  my $Q = jsc(1/(1 - $tau0), $ep);
  # s*(t) by Newton on f(s) = x tau s^2 - c s + 1, starting at s0
  my $s = jc($s0);
  for (1..6) { my $f = jadd(jsub(jsc($x, jmul($tau, jmul($s, $s))), jsc($c, $s)), jc(1));
    my $fp = jsub(jsc(2*$x, jmul($tau, $s)), jc($c)); $s = jsub($s, jmul($f, jinv($fp))); }
  my $U1 = jsub(jc(1), jmul($tau, $s));
  my $U2 = jsub(jc(($x - 1)/$x), $U1);
  ok(jeq(jadd(jsub(jsc($x, jmul($Q, jmul($U1,$U1))), jsc($x - 1, jmul($Q, $U1))), jc(-1)), jz(), $J), "LeqU at U1");
  ok(jeq(jadd(jsub(jsc($x, jmul($Q, jmul($U2,$U2))), jsc($x - 1, jmul($Q, $U2))), jc(-1)), jz(), $J), "LeqU at U2");
  # theta-polynomials P_k(y): P_1 = y, theta y^n = n (y^{n+1} - y^n); W_k = -P_k(A) - x P_k(B), A = 1/(1-U), B = 1/(1-QU)
  my $KM = 2*$M + 2; my @P; $P[1] = [R(0), R(1)];
  for my $k (2..$KM) { my $p = $P[$k-1]; my @q = map { R(0) } 0..$k;
    for my $n (1..$#$p) { next if $p->[$n]->is_zero; $q[$n+1] += $p->[$n]*$n; $q[$n] -= $p->[$n]*$n; } $P[$k] = \@q; }
  my %rows; my %Delta; my %L;
  for my $i (1, 2) { my $U = $i == 1 ? $U1 : $U2;
    my $A = jinv(jsub(jc(1), $U)); my $B = jinv(jsub(jc(1), jmul($Q, $U)));
    my (@Ap, @Bp); $Ap[0] = jc(1); $Bp[0] = jc(1);
    for my $n (1..$KM) { $Ap[$n] = jmul($Ap[$n-1], $A); $Bp[$n] = jmul($Bp[$n-1], $B); }
    my @W; for my $k (2..$KM) { my $w = jz(); for my $n (1..$k) { next if $P[$k][$n]->is_zero;
        $w = jsub($w, jsc($P[$k][$n], jadd($Ap[$n], jsc($x, $Bp[$n])))); } $W[$k] = $w; }
    # Wick (wickDirect): r^pt[m] = Sum_{n=m}^{3m} (2n)! (-1)^n/(n! 2^n W2^n) Sum_{partitions of 2n into n-m parts >= 3} prod W_k^mk/((k!)^mk mk!)
    my $iW2 = jinv($W[2]); my @iW2p = (jc(1)); push @iW2p, jmul($iW2p[-1], $iW2) for 1..3*$M;
    my @rpt;
    for my $m (0..$M) { my $tot = jz();
      if ($m == 0) { $tot = jc(1); }
      else { for my $n ($m..3*$m) { my $jn = $n - $m; next if $jn == 0;
          for my $p (parts(2*$n, $jn, 3, $KM)) { my %t; $t{$_}++ for @$p; my $mon = jc(1); my $den = R(1);
            for my $k (keys %t) { for (1..$t{$k}) { $mon = jmul($mon, $W[$k]); } $den *= fact($k)**$t{$k} * fact($t{$k}); }
            my $coef = fact(2*$n) * R((-1)**$n) / (fact($n) * R(2)**$n) / $den;
            $tot = jadd($tot, jsc($coef, jmul($mon, $iW2p[$n]))); } } }
      push @rpt, $tot; }
    # frame B: multiply by h_w(z)^-1 = exp(-Sum_k B_2k/(2k(2k-1)) z^{2k-1} (x^{1-2k} + 1 - c^{1-2k}))
    my @ex = map { R(0) } 0..$M; for my $k (1..7) { my $d = 2*$k - 1; last if $d > $M;
      $ex[$d] = -$bern[$k]/(2*$k*(2*$k-1)) * ($x**(1-2*$k) + 1 - $c**(1-2*$k)); }
    my @h = (R(1), map { R(0) } 1..$M);   # exp of ex as z-series
    for my $n (1..$M) { my $sum = R(0); $sum += $ex[$_]*$h[$n-$_]*$_ for 1..$n; $h[$n] = $sum/$n; }
    my @r0; for my $m (0..$M) { my $v = jz(); $v = jadd($v, jsc($h[$m-$_], $rpt[$_])) for 0..$m; push @r0, $v; }
    # L(U) and the code's literal Delta: Delta_1 = 1/etacan11, etacan11 = (c L2^2/x^3 - 2 L2/x^2)/(L1 - L2)^2
    my $LofU = sub { my $UU = shift; jmul(jsc(-1, jmul($Q, jadd(jsc($x, $UU), jc(1)))), jinv(jsub(jc(1), $Q))) };
    $L{$i} = $LofU->($U);
    $rows{$i} = \@r0; $rows{"pt$i"} = \@rpt; $rows{"W2_$i"} = $W[2];
  }
  my ($L1, $L2) = ($L{1}, $L{2});
  for my $i (1, 2) { my ($Li, $Lj) = $i == 1 ? ($L1, $L2) : ($L2, $L1);
    my $ec = jmul(jsub(jsc($c/$x**3, jmul($Lj, $Lj)), jsc(2/$x**2, $Lj)), jinv(jmul(jsub($Li, $Lj), jsub($Li, $Lj))));
    $Delta{$i} = jinv($ec);
    ok(jeq($Delta{$i}, jsc($x**2, $rows{"W2_$i"}), $J), "(V1) Delta_$i = x^2 W2(U_$i) (x=$x)"); }
  # H-row: r1[0] = L r0[0]; r1[m] = L r0[m] + r0[m-1]' - r0[m-1] Delta'/(2 Delta)
  my %r1; for my $i (1, 2) { my $Li = $L{$i}; my $dl = jmul(jder($Delta{$i}), jinv(jsc(2, $Delta{$i})));
    my @a = (jmul($Li, $rows{$i}[0]));
    for my $m (1..$M) { push @a, jadd(jmul($Li, $rows{$i}[$m]), jsub(jder($rows{$i}[$m-1]), jmul($rows{$i}[$m-1], $dl))); }
    $r1{$i} = \@a; }
  # (V5) second row
  for my $i (1, 2) { my ($Li, $Lj) = $i == 1 ? ($L1, $L2) : ($L2, $L1); my $dl = jmul(jder($Delta{$i}), jinv(jsc(2, $Delta{$i})));
    for my $m (1..$M) { my $lhs = jsub(jder($r1{$i}[$m-1]), jmul($r1{$i}[$m-1], $dl));
      my $rhs = jmul($Lj, jsub($r1{$i}[$m], jmul($Li, $rows{$i}[$m])));
      ok(jeq($lhs, $rhs, $J - 2), "(V5) row 2, idempotent $i, z^$m (x=$x)"); } }
  # (V6) symplecticity with the code's pairing <a,b> = x^2 (a0 b1 + a1 b0) - c x a1 b1
  my $pair = sub { my ($a0,$a1,$b0,$b1) = @_; jsub(jsc($x**2, jadd(jmul($a0,$b1), jmul($a1,$b0))), jsc($c*$x, jmul($a1,$b1))) };
  for my $i (1, 2) { for my $j (1, 2) { for my $n (0..$M) { my $s = jz();
        for my $q (0..$n) { my $p = $n - $q; my $t = $pair->($rows{$i}[$p], $r1{$i}[$p], $rows{$j}[$q], $r1{$j}[$q]);
          $s = ($q % 2) ? jsub($s, $t) : jadd($s, $t); }
        my $target = ($n == 0 && $i == $j) ? $Delta{$i} : jz();
        ok(jeq($s, $target, $J - 1), "(V6) <v_$i(w), v_$j(-w)> at w^$n (x=$x)"); } } }
  # extra: K_1 in the POINT frame too (h_w(z) h_w(-z) = 1): rebuild the H-row from the point-frame identity row
  print "point x=$x s0=$s0 tau0=$tau0 done: $ntest checks so far, $nfail failures\n";
}
print "ref_S2_jets.pl M=$M J=$J: $ntest checks, $nfail failures\n";
