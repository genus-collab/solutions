#!/usr/bin/perl
# ref_61.pl TAU HMIN HMAX [PIDX] [ENGINE_DIR]
# Referee check of Corollary 6 / Theorem 5 / Corollary 7 (mod p), independent of the author's scripts:
#  (a) Theorem 1 at n = 1 in the formal form (any rational tau):  Xi_h(j) = sigma^{h-1} j^{-2} [hb^{2h}] q_{j+1}(hb) / Pi_j,
#      q_a = Laurent coefficients of Q_j / v^j at v = oo, Q_j = exp(-j Sc(j hb d) xhat)   (j = 1..3h)
#  (b) formula (6.1): sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j-l) [lambda^{2h}] S(e1 lambda)^{1+j e2} S(e2 lambda)^{1+j e1}   (j = 0..3h)
#  (c) the engine table engine_h<h>_r1.txt if present (j = 0..3h); j = 0 of (b) vs (c) is the gap formula
#      sigma^h Int K Omega = -(2h-3)! [l^{2h}] S(e1 l) S(e2 l)  (via Xi_h(0) = (2-2h) Int K Omega, also checked).
use strict; use warnings;
my ($TAUS, $HMIN, $HMAX, $PIDX, $EDIR) = @ARGV; $PIDX //= 0;
sub isprime { my $q = shift; for (my $d = 2; $d * $d <= $q; $d++) { return 0 if $q % $d == 0; } return 1; }
my @PR; { my $q = 2147483647; while (@PR < 4) { push @PR, $q if isprime($q); $q -= 2; } }
my $P = $PR[$PIDX];
sub mm { return ($_[0] * $_[1]) % $P; }
sub ad { my $s = $_[0] + $_[1]; return $s >= $P ? $s - $P : $s; }
sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
sub pw { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub iv { my $a = shift; die "division by zero mod p" if $a % $P == 0; return pw($a, $P - 2); }
sub decmod { my $s = shift; my $neg = ($s =~ s/^-//); my $r = 0; for my $ch (split //, $s) { $r = ($r * 10 + ord($ch) - 48) % $P; } return $neg ? ng($r) : $r; }
sub ratmod { my $s = shift; $s =~ s/\s+//g; if ($s =~ m{^(-?\d+)/(\d+)$}) { return mm(decmod($1), iv(decmod($2))); } return decmod($s); }
sub im { my $k = shift; return $k >= 0 ? $k % $P : ng((-$k) % $P); }
my $tau = ratmod($TAUS); my $e2 = ad($tau, 1); my $e1 = ng($tau); my $sig = mm($tau, $e2);
my @fact = (1); push @fact, mm($fact[-1], $_) for 1..400; my @ifact = map { iv($_) } @fact;
sub binm { my ($a, $b) = @_; return 0 if $b < 0 || $b > $a; return mm($fact[$a], mm($ifact[$b], $ifact[$a-$b])); }
sub phis { my ($Nm, $K) = @_; my $half = iv(2); my (@num, @den);
  for my $i (0..$K) { $num[$i] = mm(pw(mm($Nm, $half), 2*$i+1), $ifact[2*$i+1]); $den[$i] = mm(pw($half, 2*$i+1), $ifact[2*$i+1]); }
  my @q; my $id0 = iv($den[0]); for my $i (0..$K) { my $s = $num[$i]; for my $j (1..$i) { $s = ad($s, ng(mm($den[$j], $q[$i-$j]))); } $q[$i] = mm($s, $id0); } return \@q; }
# q_{mu+1} coefficient of hb^{2h} only needed: compute Laurent series up to a = mu+1 with hb-degree <= 2h
sub qcoef { my ($mu, $HTq) = @_; my $A = $mu + 1; my $K = int($HTq / 2);
  my $Na = mm($tau, im($mu)); my $Nb = mm($e2, im($mu)); my $pa = phis($Na, $K); my $pb = phis($Nb, $K);
  my @L = map { [ (0) x ($HTq + 1) ] } 0..$A;
  for my $a (1..$A) { my $c = mm($Nb, iv($a)); $c = ng($c) if $a % 2 == 0; $L[$a][0] = ad($L[$a][0], $c); }
  for my $k (1..$K) { my $f = $fact[2*$k-1]; my $ca = mm($f, mm($pa->[$k], iv(pw($tau, 2*$k)))); my $cb = mm($f, mm($pb->[$k], iv(pw($e2, 2*$k))));
    $L[2*$k][2*$k] = ad($L[2*$k][2*$k], $ca) if 2*$k <= $A;
    for my $b (0..$A - 2*$k) { my $cf = binm(2*$k + $b - 1, $b); $cf = ng($cf) if $b % 2 == 1; $L[2*$k + $b][2*$k] = ad($L[2*$k + $b][2*$k], ng(mm($cb, $cf))); } }
  my @E = ([ 1, (0) x $HTq ]);
  for my $a (1..$A) { my @acc = (0) x ($HTq + 1);
    for my $j (1..$a) { my $Lj = $L[$j]; my $Ea = $E[$a-$j]; for my $p (0..$HTq) { next unless $Lj->[$p]; my $c = mm($Lj->[$p], $j);
        for my $q (0..$HTq - $p) { next unless $Ea->[$q]; $acc[$p+$q] = ($acc[$p+$q] + $c * $Ea->[$q]) % $P; } } }
    my $ia = iv($a); $E[$a] = [ map { mm($_, $ia) } @acc ]; }
  return $E[$A]; }
sub Pi { my $mu = shift; my $r = 1; for my $a (1..$mu-1) { $r = mm($r, ad(mm(im($mu), $tau), im($a))); } return mm($r, $ifact[$mu-1]); }
# log S(y) series in y^{2k}: S = 1/Sc, Sc(y) = sum (y/2)^{2i}/(2i+1)!
my $KMAX = $HMAX + 1;
my @sc; for my $i (0..$KMAX) { $sc[$i] = mm(pw(iv(2), 2*$i), $ifact[2*$i+1]); }   # coefficients in Y = y^2
my @lsc = (0); for my $n (1..$KMAX) { my $s = mm($n, $sc[$n]); for my $k (1..$n-1) { $s = ad($s, ng(mm(mm($k, $lsc[$k]), $sc[$n-$k]))); } $lsc[$n] = mm($s, iv($n)); }
my @lS = map { ng($_) } @lsc;   # log S in powers of y^2
sub f61 { my ($j, $h) = @_;   # sigma^h Xi_h(j) by (6.1)
  my $a = ad(1, mm(im($j), $e2)); my $b = ad(1, mm(im($j), $e1));
  my @ex = (0); for my $k (1..$h) { $ex[$k] = ad(mm($a, mm($lS[$k], pw($e1, 2*$k))), mm($b, mm($lS[$k], pw($e2, 2*$k)))); }
  my @r = (1); for my $n (1..$h) { my $s = 0; for my $m (1..$n) { $s = ad($s, mm(mm($m, $ex[$m]), $r[$n-$m])); } $r[$n] = mm($s, iv($n)); }
  my $pre = 1; $pre = mm($pre, im($j - $_)) for 1..(2*$h-2); return mm($pre, $r[$h]); }
sub gapf { my $h = shift;  # sigma^h Int K Omega by the gap formula
  my @ex = (0); for my $k (1..$h) { $ex[$k] = ad(mm($lS[$k], pw($e1, 2*$k)), mm($lS[$k], pw($e2, 2*$k))); }
  my @r = (1); for my $n (1..$h) { my $s = 0; for my $m (1..$n) { $s = ad($s, mm(mm($m, $ex[$m]), $r[$n-$m])); } $r[$n] = mm($s, iv($n)); }
  return ng(mm($fact[2*$h-3], $r[$h])); }
sub load_engine { my ($file, $h) = @_; my %E; open my $fh, '<', $file or return undef;
  while (my $line = <$fh>) { chomp $line; my @f = split /\s*\|\s*/, $line; next unless @f >= 4; next unless $f[0] == $h && $f[1] == 1;
    for my $t (split /\s*;\s*/, $f[3]) { next unless $t =~ /^(\d+):(\S+)$/; push @{$E{$f[2]}}, [ $1, ratmod($2) ]; } } close $fh; return \%E; }
sub xi_eng { my ($E, $j) = @_; my $mis = ng(iv($sig)); my $v = 0;
  for my $k (keys %$E) { my $s = 0; for my $t (@{$E->{$k}}) { $s = ad($s, mm($t->[1], pw(im($j), $t->[0]))); } $v = ad($v, mm($s, pw($mis, $k))); } return $v; }
my ($nab, $bab, $nbc, $bbc, $ngap, $bgap) = (0, 0, 0, 0, 0, 0);
for my $h ($HMIN..$HMAX) {
  my $E; if (defined $EDIR) { my $f = ($h == 7) ? "$EDIR/../round2/R2_onepoint/engine_h7_r1.txt" : "$EDIR/engine_h${h}_r1.txt"; $E = load_engine($f, $h) if -e $f; }
  my ($la, $lb, $lc) = (0, 0, 0);
  for my $j (0..3*$h) {
    my $b = mm(f61($j, $h), iv(pw($sig, $h)));   # Xi_h(j) by (6.1)
    if ($j >= 1 && Pi($j) == 0) { my $q = qcoef($j, 2*$h); if ($q->[2*$h] != 0) { $bab++; print "  INCONSISTENT: Pi_j = 0 but [x^j]F != 0 (tau=$TAUS h=$h j=$j)\n"; } }
    if ($j >= 1 && Pi($j) != 0) { my $q = qcoef($j, 2*$h); my $F = mm($q->[2*$h], iv(mm($j, $j))); my $a = mm(pw($sig, $h - 1), mm($F, iv(Pi($j))));
      $nab++; $la++; if ($a != $b) { $bab++; print "  MISMATCH Thm1 vs (6.1): tau=$TAUS h=$h j=$j\n"; } }
    if ($E) { my $c = xi_eng($E, $j); $nbc++; $lc++; if ($c != $b) { $bbc++; print "  MISMATCH (6.1) vs engine: tau=$TAUS h=$h j=$j\n"; } }
  }
  if ($h >= 2) { my $gap = mm(gapf($h), iv(pw($sig, $h)));   # Int K Omega
    my $x0 = mm(f61(0, $h), iv(pw($sig, $h)));                # Xi_h(0) by (6.1)
    $ngap++; if (mm(im(2 - 2*$h), $gap) != $x0) { $bgap++; print "  MISMATCH gap: tau=$TAUS h=$h\n"; } }
  print "tau=$TAUS h=$h: Thm1(n=1) vs (6.1) at $la values" . ($E ? ", (6.1) vs engine at $lc values" : ", no engine table") . "\n"; }
print "SUMMARY tau=$TAUS h=$HMIN..$HMAX p=$P: Thm1 vs (6.1): $nab compared, $bab mismatches; (6.1) vs engine: $nbc compared, $bbc mismatches; gap vs (6.1) at j=0: $ngap compared, $bgap mismatches\n";
