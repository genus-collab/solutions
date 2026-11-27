#!/usr/bin/perl
# Referee check for S2.4(d)/(I1): the code's frame-B singular row, built as in rc_lib.m (Wick series in y = log U at U_1,
# W_k = (U d/dU)^k W, twist h_w^-1, code's literal Delta_1, H-row recursion), as exact t-series at rational x, against
#   KEY FORMULA  G = rho r_{0,1} = c/(c+z) (tau/t) 2F1(1/z+1, 1/z+1; c/z+2; tau) Exp[-dphi/z],  rho = (Delta_1 t^2/(-x c^3))^(-1/2)
#   (I1)         t rho r_{1,1} = (t L_1 + z) G + z t dG/dt
# Trick: t^k W_k is a power series in t, so every Wick monomial (t^2 W_2)^-n prod (t^k W_k) is a power series (sum k = 2n).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
my $M = @ARGV ? $ARGV[0] : 5;      # z-order
my $K = @ARGV > 1 ? $ARGV[1] : 4;  # t-order compared
my $N = $K + 2;                     # internal t-length
my ($ntest, $nfail) = (0, 0);
sub ok { my ($c, $msg) = @_; $ntest++; unless ($c) { $nfail++; print "FAIL: $msg\n"; } }
sub z0 { [ map { R(0) } 0..$N ] }
sub c0 { my $a = z0(); $a->[0] = R($_[0]); $a }
sub add { my ($a,$b)=@_; [ map { $a->[$_] + $b->[$_] } 0..$N ] }
sub sub_ { my ($a,$b)=@_; [ map { $a->[$_] - $b->[$_] } 0..$N ] }
sub sc { my ($s,$a)=@_; [ map { $a->[$_]*$s } 0..$N ] }
sub mul { my ($a,$b)=@_; my $c = z0(); for my $i (0..$N) { next if $a->[$i]->is_zero; for my $k (0..$N-$i) { $c->[$i+$k] += $a->[$i]*$b->[$k]; } } $c }
sub inv { my $a = shift; die "inv0" if $a->[0]->is_zero; my $b = z0(); $b->[0] = 1/$a->[0]; for my $n (1..$N) { my $s = R(0); $s += $a->[$_]*$b->[$n-$_] for 1..$n; $b->[$n] = -$s/$a->[0]; } $b }
sub sqrt1 { my $a = shift; die "sqrt1" unless $a->[0] == 1; my $b = z0(); $b->[0] = R(1); for my $n (1..$N) { my $s = R(0); $s += $b->[$_]*$b->[$n-$_] for 1..$n-1; $b->[$n] = ($a->[$n]-$s)/2; } $b }
sub der { my $a = shift; my $b = z0(); $b->[$_-1] = $a->[$_]*$_ for 1..$N; $b }
sub integ { my $a = shift; my $b = z0(); $b->[$_] = $a->[$_-1]/$_ for 1..$N; $b }
sub up { my $a = shift; my $b = z0(); $b->[$_] = $a->[$_-1] for 1..$N; $b }
sub tpow { my $k = shift; my $a = z0(); $a->[$k] = R(1) if $k <= $N; $a }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
sub parts { my ($S, $k, $lo, $hi) = @_; my @out; return ([]) if $k == 0 && $S == 0; return () if $k == 0;
  for my $p (reverse $lo..($hi < $S ? $hi : $S)) { next if $p*$k < $S; next if $lo*($k-1) > $S - $p;
    for my $rest (parts($S - $p, $k - 1, $lo, $p)) { push @out, [$p, @$rest]; } } @out }
my @bern = (undef, R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"));

for my $x (R(2), R("3/5"), R(-3)) { my $c = 1 + $x;
  my $tt = z0(); $tt->[$_] = R((($_ % 2) ? -1 : 1))/fact($_ + 1) for 0..$N;         # tau/t
  my $tau = up($tt);
  my $Q = z0(); $Q->[$_] = R(1)/fact($_) for 0..$N;                                    # e^t
  my $qm1 = z0(); $qm1->[$_] = R(1)/fact($_ + 1) for 0..$N;                           # (Q - 1)/t
  my $dc = sqrt1(sub_(c0(1), sc(4*$x/$c**2, $tau)));                                   # delta/c (principal branch)
  my $ss = sc(1/$c, sc(R(2), inv(add(c0(1), $dc))));                                   # s* = 2/(c + delta)
  my $U1 = sub_(c0(1), mul($tau, $ss));                                                # code's U1 = (S0 + sq)/2
  # a = t/(1 - U1) = 1/(tt s*),  b = t/(1 - Q U1) = -1/(Q tt (1 - s*))
  my $a = inv(mul($tt, $ss)); my $b = sc(-1, inv(mul($Q, mul($tt, sub_(c0(1), $ss)))));
  my $KM = 2*$M + 2; my @P; $P[1] = [R(0), R(1)];
  for my $k (2..$KM) { my $p = $P[$k-1]; my @q = map { R(0) } 0..$k;
    for my $n (1..$#$p) { next if $p->[$n]->is_zero; $q[$n+1] += $p->[$n]*$n; $q[$n] -= $p->[$n]*$n; } $P[$k] = \@q; }
  my (@ap, @bp); $ap[0] = c0(1); $bp[0] = c0(1);
  for my $n (1..$KM) { $ap[$n] = mul($ap[$n-1], $a); $bp[$n] = mul($bp[$n-1], $b); }
  my @Wt; for my $k (2..$KM) { my $w = z0();                                           # Wt_k = t^k W_k(U1)
    for my $n (1..$k) { next if $P[$k][$n]->is_zero; $w = sub_($w, sc($P[$k][$n], mul(tpow($k - $n), add($ap[$n], sc($x, $bp[$n]))))); }
    $Wt[$k] = $w; }
  my $iW2 = inv($Wt[2]); my @iW2p = (c0(1)); push @iW2p, mul($iW2p[-1], $iW2) for 1..3*$M;
  my @rpt = (c0(1));
  for my $m (1..$M) { my $tot = z0();
    for my $n ($m..3*$m) { my $jn = $n - $m; next if $jn == 0;
      for my $p (parts(2*$n, $jn, 3, $KM)) { my %t; $t{$_}++ for @$p; my $mon = c0(1); my $den = R(1);
        for my $k (keys %t) { for (1..$t{$k}) { $mon = mul($mon, $Wt[$k]); } $den *= fact($k)**$t{$k} * fact($t{$k}); }
        $tot = add($tot, sc(fact(2*$n)*R((-1)**$n)/(fact($n)*R(2)**$n)/$den, mul($mon, $iW2p[$n]))); } }
    push @rpt, $tot; }
  # frame B
  my @ex = map { R(0) } 0..$M; for my $k (1..7) { my $d = 2*$k - 1; last if $d > $M;
    $ex[$d] = -$bern[$k]/(2*$k*(2*$k-1)) * ($x**(1-2*$k) + 1 - $c**(1-2*$k)); }
  my @h = (R(1), map { R(0) } 1..$M); for my $n (1..$M) { my $s = R(0); $s += $ex[$_]*$h[$n-$_]*$_ for 1..$n; $h[$n] = $s/$n; }
  my @r0; for my $m (0..$M) { my $v = z0(); $v = add($v, sc($h[$m-$_], $rpt[$_])) for 0..$m; push @r0, $v; }
  # L1, L2, code's literal Delta_1 (times t^2), rho
  my $tL1 = mul(mul($Q, add(sc($x, $U1), c0(1))), inv($qm1));                         # t L1 = -Q (x U1 + 1) t/(1 - Q)
  my $L2 = mul(mul($Q, sc($x, mul($tt, $ss))), inv($qm1));                            # L2 = -Q (x U2 + 1)/(1 - Q), x U2 + 1 = x tau s*
  my $d12 = sub_($tL1, up($L2));                                                       # t (L1 - L2)
  my $tD1 = mul(sc($x**3, mul($d12, $d12)), inv(mul($L2, sub_(sc($c, $L2), c0(2*$x)))));   # t^2 Delta_1
  my $rho = inv(sqrt1(sc(1/(-$x*$c**3), $tD1)));
  # closed form, bivariate: CF[k][m + K] = [t^k][z^m], m in [-K, M]
  my $ZL = $M + $K;                                                                    # z-length of power series needed
  my $zmul = sub { my ($p, $q) = @_; my @r = map { R(0) } 0..$ZL; for my $i (0..$ZL) { next if $p->[$i]->is_zero; for my $j (0..$ZL-$i) { $r[$i+$j] += $p->[$i]*$q->[$j]; } } \@r };
  my $zinvlin = sub { my ($A0, $A1) = @_; [ map { (1/$A0) * (-$A1/$A0)**$_ } 0..$ZL ] };   # 1/(A0 + A1 z)
  my $dphi = integ(sub_(c0(1), sc($x, $ss))); my $dpt = [ @{$dphi}[1..$N], R(0) ];       # dphi/t
  my @CF = map { [ map { R(0) } 0..($K + $M) ] } 0..$K;
  my $cz = $zinvlin->($c, R(1)); $cz = [ map { $_*$c } @$cz ];                           # c/(c+z)
  for my $n (0..$K) { my $Fh = [ R(1), map { R(0) } 1..$ZL ];
    for my $j (1..$n) { $Fh = $zmul->($Fh, $zmul->([R(1), R($j), map { R(0) } 2..$ZL], [R(1), R($j), map { R(0) } 2..$ZL])); }
    for my $j (2..$n+1) { $Fh = $zmul->($Fh, $zinvlin->($c, R($j))); }
    $Fh = $zmul->($Fh, $cz); $Fh = [ map { $_/fact($n) } @$Fh ];
    my $ttn = c0(1); $ttn = mul($ttn, $tt) for 1..$n;
    my $dpj = c0(1);
    for my $jj (0..$K - $n) {
      my $tser = mul(mul($tt, $ttn), sc(R((-1)**$jj)/fact($jj), $dpj));               # (tau/t) tt^n (-dphi/t)^jj/jj!
      for my $kk (0..$K - $n - $jj) { my $coef = $tser->[$kk]; next if $coef->is_zero; my $k = $kk + $n + $jj;
        for my $i (0..$ZL) { my $m = $i - $n - $jj; next if $m > $M; $CF[$k][$m + $K] += $coef*$Fh->[$i]; } }
      $dpj = mul($dpj, $dpt); } }
  for my $k (0..$K) {
    for my $m (-$K..-1) { ok($CF[$k][$m + $K]->is_zero, "closed form: no z^$m at t^$k (x=$x)"); }
    for my $m (0..$M) { my $Gk = mul($rho, $r0[$m])->[$k];
      ok($Gk == $CF[$k][$m + $K], "KEY FORMULA [t^$k z^$m] (x=$x): code $Gk vs closed form $CF[$k][$m+$K]"); } }
  # (I1): t rho r11 = (t L1 + z) G + z t G', with r11 from the code's H-row (H)
  my $dl = sub_(mul(der($tD1), inv($tD1)), c0(0));                                    # (t^2 D1)'/(t^2 D1) = D1'/D1 + 2/t
  # t D1'/(2 D1) = t (t^2D1)'/(2 t^2 D1) - 1
  my $tdl = sub_(sc(R("1/2"), up($dl)), c0(1));
  my @tr1; for my $m (0..$M) { my $v = mul($tL1, $r0[$m]);
    if ($m >= 1) { $v = add($v, sub_(up(der($r0[$m-1])), mul($r0[$m-1], $tdl))); } push @tr1, $v; }
  my @G = map { mul($rho, $r0[$_]) } 0..$M;
  for my $m (0..$M) { my $lhs = mul($rho, $tr1[$m]);
    my $rhs = mul($tL1, $G[$m]); if ($m >= 1) { $rhs = add($rhs, add($G[$m-1], up(der($G[$m-1])))); }
    for my $k (0..$K-1) { ok($lhs->[$k] == $rhs->[$k], "(I1) [t^$k z^$m] (x=$x)"); } }
  print "x=$x: rho = ", join(", ", map { $rho->[$_] } 0..2), "...; $ntest checks so far, $nfail failures\n";
}
print "ref_S2_key.pl M=$M K=$K: $ntest checks, $nfail failures\n";
