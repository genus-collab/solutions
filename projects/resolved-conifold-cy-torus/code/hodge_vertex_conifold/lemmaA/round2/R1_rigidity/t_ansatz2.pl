#!/usr/bin/perl
# t_ansatz2.pl -- test of the deformed-curve formula for the genus-G linearized solution at a general point of the slice:
#   F^Y(P) = [lambda^{-1}] Psi(lambda) log( z_P(lambda) u(lambda) ),   Psi(lambda) = sum_c Y(c) u(lambda)^c,  u = (lambda-1)^tau lambda^{-1-tau},
#   z_P(lambda) solves  Y_P(z) - sum_{j>=2} P_j z^j = lambda + tau,   Y_P(z) = W_* Y_0(z/W_*),  W_* = 1 + sum_j Pi_j P_j W_*^j.
# Along the line P_j = eps c_j (c_j random integers, 2 <= j <= J), the eps-expansion is compared with the KP propagation
# (kp_dump.pl: multi-leaf values Y(T) as linear forms in the one-leaf values) for random one-leaf values.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M, $J, $seed) = @ARGV; $taus //= '2/7'; $M //= 3; $J //= 4; $seed //= 5; my $tau = R($taus); my $sig = $tau*(1+$tau);
srand($seed); my %cj; $cj{$_} = R(int(rand(7)) - 3) for 2..$J;
my $Lmax = $M*($J-1) + 1;   # highest level reached by eps^M terms
my @yv = (R(0)); push @yv, R(int(rand(11)) - 5) / R(1 + int(rand(3))) for 1..$Lmax+2;   # y_L = Y(L+1), y_0 = 0
$yv[1] = R(0);   # Xi(2) = 0 (Theorem A) -- not needed for the test but harmless
# ---------------- ansatz side
my $NEG = $Lmax + 3; my $OFF = $NEG; my $TOPP = $M*$J + 3; my $LEN = $OFF + $TOPP;
sub lz { [ map { R(0) } 0..$LEN ] }
sub lmul { my ($a, $b) = @_; my $c = lz(); for my $i (0..$LEN) { next if $a->[$i]->is_zero; for my $j (0..$LEN) { next if $b->[$j]->is_zero; my $k = $i + $j - $OFF; next if $k < 0 || $k > $LEN; $c->[$k] = $c->[$k] + $a->[$i]*$b->[$j]; } } return $c; }
sub ladd { my ($a, $b, $s) = @_; $s //= R(1); return [ map { $a->[$_] + $s*$b->[$_] } 0..$LEN ]; }
sub lmon { my ($p, $co) = @_; my $v = lz(); my $i = $p + $OFF; $v->[$i] = $co // R(1) if $i >= 0 && $i <= $LEN; return $v; }
sub bz { [ map { lz() } 0..$M ] }
sub bmul { my ($a, $b) = @_; my $c = bz(); for my $i (0..$M) { for my $j (0..$M-$i) { $c->[$i+$j] = ladd($c->[$i+$j], lmul($a->[$i], $b->[$j])); } } return $c; }
sub badd { my ($a, $b, $s) = @_; return [ map { ladd($a->[$_], $b->[$_], $s) } 0..$M ]; }
sub bconst { my $v = bz(); $v->[0] = lmon(0, $_[0]); return $v; }
# W_*(eps)
my @W = (R(1), map { R(0) } 1..$M);
for my $it (1..$M+1) { my @nw = (R(1), map { R(0) } 1..$M); my $wp = [@W];
  for my $j (2..$J) { $wp = series_mul($wp, \@W, $M); my @Pjs = (R(0), $cj{$j}, map { R(0) } 2..$M); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @Pjs ], $wp, $M); $nw[$_] = $nw[$_] + $q->[$_] for 0..$M; } @W = @nw; }
my $z0 = lz(); for my $i (0..$OFF+1) { my $b = R(1); for my $q (0..$i-1) { $b = $b * ($tau + $q) / ($q + 1); } $z0 = ladd($z0, lmon(1-$i, $b)); }
my $alB = bz(); for my $m (0..$M) { $alB->[$m] = lmon(0, $W[$m]); }
sub binv_lead { my $zz = shift; my $q = bz(); for my $m (0..$M) { $q->[$m] = lmul($zz->[$m], lmon(-1)); } my $d = badd($q, bconst(R(1)), R(-1));
  my $inv = bconst(R(1)); my $pw = bconst(R(1)); for my $k (1..$OFF + $TOPP + 2) { $pw = bmul($pw, $d); last unless grep { grep { !$_->is_zero } @$_ } @$pw; $inv = badd($inv, $pw, R((-1)**$k)); }
  my $r = bz(); for my $m (0..$M) { $r->[$m] = lmul($inv->[$m], lmon(-1)); } return $r; }
my $zi00 = binv_lead([ $z0, map { lz() } 1..$M ])->[0];
my $dY = lmon(0, R(1)); { my $p = lmon(0, R(1)); for my $k (1..$OFF + 2) { $p = lmul($p, $zi00); $dY = ladd($dY, lmul($p, $zi00), $sig*Pi_at($k,$tau)*$k/($k+1)); } }
my $idY; { my $dd = ladd($dY, lmon(0, R(1)), R(-1)); $idY = lmon(0, R(1)); my $pp = lmon(0, R(1)); for my $k (1..$OFF+$TOPP) { $pp = lmul($pp, $dd); last unless grep { !$_->is_zero } @$pp; $idY = ladd($idY, $pp, R((-1)**$k)); } }
my $z = bz(); $z->[0] = [ @$z0 ];
for my $it (1..2*$M+3) {
  my $zi = binv_lead($z); my $Fz = [ map { [ @$_ ] } @$z ]; my $zk = bconst(R(1)); my $alp = [ map { [ @$_ ] } @$alB ];
  for my $k (1..$OFF + 2) { $zk = bmul($zk, $zi); $alp = bmul($alp, $alB); $Fz = badd($Fz, bmul($alp, $zk), -$sig*Pi_at($k,$tau)/($k+1)); }
  my $zj = [ map { [ @$_ ] } @$z ]; for my $j (2..$J) { $zj = bmul($zj, $z); my $Pb = bz(); $Pb->[1] = lmon(0, $cj{$j}); $Fz = badd($Fz, bmul($Pb, $zj), R(-1)); }
  $Fz = badd($Fz, bconst($tau), R(-1)); $Fz->[0] = ladd($Fz->[0], lmon(1, R(-1)));
  my $corr = bz(); for my $m (0..$M) { $corr->[$m] = lmul($Fz->[$m], $idY); } $z = badd($z, $corr, R(-1)); }
my $ratio = bz(); for my $m (0..$M) { $ratio->[$m] = lmul($z->[$m], $zi00); }
my $sm = badd($ratio, bconst(R(1)), R(-1)); my $K = bz(); { my $pw = bconst(R(1)); for my $k (1..$M) { $pw = bmul($pw, $sm); $K = badd($K, $pw, R((-1)**($k+1))/$k); } }
# Psi(lambda) = sum_c y_{c-1} u^c, u = 1/z_0
my $Psi = lz(); { my $up = lmon(0, R(1)); for my $cc (1..$Lmax+2) { $up = lmul($up, $zi00); $Psi = ladd($Psi, $up, $yv[$cc-1]); } }
my @ans; for my $m (1..$M) { my $pr = lmul($Psi, $K->[$m]); push @ans, $pr->[-1 + $OFF]; }
# ---------------- KP side
my $out = `perl "$FindBin::Bin/kp_dump.pl" 1 $Lmax $taus`; my %Yf; for (split /\n/, $out) { next unless /^Y (\S+) : (.*)$/; $Yf{$1} = [ map { R($_) } split ' ', $2 ]; }
sub Yv { my @T = @_; my $f = $Yf{join(',', sort { $b <=> $a } @T)}; return R(0) unless $f; my $v = R(0); for my $i (0..$#$f) { $v = $v + $f->[$i]*$yv[$i+1]; } return $v; }
my @kp = map { R(0) } 1..$M;
sub msets { my ($r, $min) = @_; return ([]) if $r == 0; my @o; for my $j ($min..$J) { push @o, [$j, @$_] for msets($r-1, $j); } return @o; }
for my $r (1..$M) { for my $T (msets($r, 2)) { my $w = R(1); $w = $w * $cj{$_} for @$T; next if $w->is_zero; my $lv = 0; $lv += $_ - 1 for @$T; next if $lv > $Lmax;
    $kp[$r-1] = $kp[$r-1] + $w * Yv(@$T) / aut(@$T); } }
my $bad = 0; for my $m (1..$M) { my $ok = $ans[$m-1] == $kp[$m-1]; $bad++ unless $ok; print "eps^$m: ansatz $ans[$m-1]   KP $kp[$m-1]   ", ($ok ? "OK" : "DIFF"), "\n"; }
print "t_ansatz2 (tau=$taus, M=$M, J=$J): ", ($bad ? "$bad orders differ" : "agreement"), "\n";
