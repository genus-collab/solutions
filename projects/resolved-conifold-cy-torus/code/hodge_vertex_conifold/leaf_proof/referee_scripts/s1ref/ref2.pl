#!/usr/bin/perl
# Referee check R5 for S1 s.7.4: transport recursion P Phi_m' + R Phi_m + T Phi_{m-1} = 0 (Phi_0(0) = 1, Phi_m(0) = 0)
# versus [z^m] (Exp[-dphi/z] F) computed directly as tau-series (Laurent bookkeeping), m = 0..MM, mod tau^K.
# Also R6: homogeneity -- w_k := tau^k W_k(U1) is analytic at tau = 0 with w_k(0) = (-1)^k phi_0^(k)(1/c) (k = 2..12).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
our $K;
sub z0 { [ map { R(0) } 1..$K ] }
sub cst { my $a = shift; my $s = z0(); $s->[0] = ref $a ? $a->copy : R($a); $s }
sub add { my ($a, $b) = @_; [ map { $a->[$_] + $b->[$_] } 0..$K-1 ] }
sub sbt { my ($a, $b) = @_; [ map { $a->[$_] - $b->[$_] } 0..$K-1 ] }
sub scl { my ($c, $a) = @_; [ map { $c * $a->[$_] } 0..$K-1 ] }
sub mul { my ($a, $b) = @_; my $r = z0();
  for my $i (0..$K-1) { next if $a->[$i]->is_zero; for my $j (0..$K-1-$i) { $r->[$i+$j] += $a->[$i] * $b->[$j] } } $r }
sub inv { my $a = shift; die "inv" if $a->[0]->is_zero; my $r = z0(); $r->[0] = 1 / $a->[0];
  for my $n (1..$K-1) { my $s = R(0); $s += $a->[$_] * $r->[$n-$_] for 1..$n; $r->[$n] = -$s / $a->[0] } $r }
sub sqrt1 { my $a = shift; die "sqrt1" unless $a->[0] == 1; my $r = z0(); $r->[0] = R(1);
  for my $n (1..$K-1) { my $s = $a->[$n]->copy; $s -= $r->[$_] * $r->[$n-$_] for 1..$n-1; $r->[$n] = $s / 2 } $r }
sub integ { my $a = shift; my $r = z0(); $r->[$_] = $a->[$_-1] / $_ for 1..$K-1; $r }
sub der { my $a = shift; my $r = z0(); $r->[$_-1] = $a->[$_] * $_ for 1..$K-1; $r }   # top coefficient lost (becomes 0)
sub pw { my ($a, $n) = @_; my $r = cst(1); $r = mul($r, $a) for 1..$n; $r }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
my $bad = 0; sub chk { my ($ok, $msg) = @_; print(($ok ? "ok    " : "FAIL  "), $msg, "\n"); $bad++ unless $ok }
for my $x (R(2), R("1/3"), R("-5/2")) {
  my $c = 1 + $x; my $MM = 3; my $KT = 7; $K = $KT + 2 * $MM + 2;   # extra tau-orders: each transport step loses accuracy via T (2 derivatives)
  my $T = z0(); $T->[1] = R(1);
  my $dc = sqrt1(add(cst(1), scl(-4 * $x / $c**2, $T)));
  my $delta = scl($c, $dc);
  my $sst = scl(2 / $c, inv(add(cst(1), $dc)));
  my $alg = mul(sqrt1($dc), scl(R(1)/2, add(cst(1), $dc)));
  my $P = $delta;
  my $Rr = scl(-$x, mul(add(cst($c), scl(3, $delta)), inv(mul($delta, add(cst($c), $delta)))));
  my $Top = sub { my $f = shift; my $f1 = der($f); my $f2 = der($f1);
    sbt(add(mul(mul($T, sbt(cst(1), $T)), $f2), mul(sbt(cst(2), scl(3, $T)), $f1)), $f) };
  my @Phi; $Phi[0] = inv($alg);
  # check transport m = 0: P Phi_0' + R Phi_0 = 0
  my $t0 = add(mul($P, der($Phi[0])), mul($Rr, $Phi[0]));
  my $ok0 = 1; for my $i (0..$K-2) { $ok0 = 0 unless $t0->[$i]->is_zero } chk($ok0, "x=$x: Phi_0 = 1/alg solves P Phi_0' + R Phi_0 = 0");
  for my $m (1..$MM) { my $rhs = mul($Top->($Phi[$m-1]), inv(mul($P, $Phi[0])));
    $Phi[$m] = scl(-1, mul($Phi[0], integ($rhs))) }
  # direct: Exp[-dphi/z] F, [z^m] for m = 0..MM, mod tau^KT
  my $Kd = $KT; local $K = $Kd;
  my $Ts = z0(); $Ts->[1] = R(1);
  my $dcs = sqrt1(add(cst(1), scl(-4 * $x / $c**2, $Ts)));
  my $ss = scl(2 / $c, inv(add(cst(1), $dcs)));
  my $dphi = integ(mul($ss, inv(sbt(cst(1), mul($Ts, $ss)))));
  my $ZT = $MM + $Kd + 1; my %F;
  for my $n (0..$Kd-1) { my @S = (R(1)); push @S, R(0) for 1..$ZT+$n;
    my $pmul = sub { my ($q) = @_; my @r = map { R(0) } 0..$#S; for my $i (0..$#S) { next if $S[$i]->is_zero; for my $j (0..$#$q) { last if $i+$j > $#S; $r[$i+$j] += $S[$i] * $q->[$j] } } @S = @r };
    for my $j (1..$n) { $pmul->([R(1), R($j)]); $pmul->([R(1), R($j)]) }
    for my $j (2..$n+1) { $pmul->([ map { (R(-$j)/$c)**$_ / $c } 0..$#S ]) }
    for my $i (0..$#S) { my $zp = $i - $n; next if $zp > $ZT; my $ts = z0(); $ts->[$n] = $S[$i] / fact($n);
      $F{$zp} = defined $F{$zp} ? add($F{$zp}, $ts) : $ts } }
  my %E; { my $q = cst(1); for my $l (0..$Kd-1) { $E{-$l} = scl(R(1)/fact($l), $q); $q = mul($q, scl(-1, $dphi)) } }
  my %Z; for my $p (keys %F) { for my $q (keys %E) { my $r = $p + $q; next if $r > $MM; my $pr = mul($F{$p}, $E{$q}); $Z{$r} = defined $Z{$r} ? add($Z{$r}, $pr) : $pr } }
  my $okneg = 1; for my $p (keys %Z) { next unless $p < 0; for (@{$Z{$p}}) { $okneg = 0 unless $_->is_zero } }
  chk($okneg, "x=$x: Exp[-dphi/z] F has no negative z-powers mod tau^$Kd");
  for my $m (0..$MM) { my $ok = 1; for my $i (0..$Kd-1) { $ok = 0 unless $Z{$m}[$i] == $Phi[$m][$i] }
    chk($ok, "x=$x: [z^$m] (Exp[-dphi/z] F) = Phi_$m of the transport recursion, mod tau^$Kd") }
}
# R6 homogeneity
for my $x (R(2), R("1/3"), R("-5/2")) { my $c = 1 + $x; local $K = 4;
  my $T = z0(); $T->[1] = R(1);
  my $dc = sqrt1(add(cst(1), scl(-4 * $x / $c**2, $T)));
  my $sst = scl(2 / $c, inv(add(cst(1), $dc)));
  my $a = scl(-1, mul(sbt(cst(1), $T), inv(sbt(cst(1), $sst)))); my $b = inv($sst);
  my @P = ([R(0), R(1)]);
  for my $n (1..12) { my @d = map { $P[$n-1][$_] * $_ } 1..$#{$P[$n-1]}; my @r = map { R(0) } 0..$#d + 2;
    for my $i (0..$#d) { $r[$i+2] += $d[$i]; $r[$i+1] -= $d[$i] } push @P, \@r }
  my $ok = 1;
  for my $k (2..12) { my $s = z0();
    for my $i (0..$#{$P[$k-1]}) { my $co = $P[$k-1][$i]; next if $co->is_zero; my $tp = pw($T, $k - $i);
      $s = sbt($s, scl($x * $co, mul(pw($a, $i), $tp))); $s = sbt($s, scl($co, mul(pw($b, $i), $tp))) }
    my $phik = fact($k-1) * $c**$k * (R(-1)**($k-1) - $x**(1-$k));
    $ok = 0 unless $s->[0] == R(-1)**$k * $phik }
  chk($ok, "x=$x: tau^k W_k(U1) -> (-1)^k phi_0^(k)(1/c) for k = 2..12 (each Wick monomial is analytic at t = 0)") }
print "\n", ($bad ? "$bad FAILURES" : "all referee checks passed"), "\n";
