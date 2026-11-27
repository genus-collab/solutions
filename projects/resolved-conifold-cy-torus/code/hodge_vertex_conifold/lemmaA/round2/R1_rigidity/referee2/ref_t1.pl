#!/usr/bin/perl
# ref_t1.pl -- referee 2: independent test of T1 (closed form of the polar-order relations) and of the internal
# consistency of the KP propagation.
#   relation_k := [t^k] sum_T Y(T) P(t)^T/|Aut T|  on the leaf curve P_j(t) = (-1)^j (Pi_j/j) t^(j-1) H_j/H_1^j,
#   T1: relation_k = kappa_k sum_{L<=k} y_L (-1)^(k-L) [(1+tau) binom((L+1)tau, k-L) - tau binom((L+1)(1+tau)+k-L-1, k-L)].
# Usage: perl ref_t1.pl K tau mode [strategy] [prime]
#   mode = sym  : one-leaf values are formal unknowns (coefficient-wise comparison)
#   mode = rN   : N random one-leaf vectors (seeded)
#   strategy    : big | small | mix (pair used in Theorem D for each T); 'cons' compares big vs small vs mix
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rk;
my ($K, $ts, $mode, $strat, $prime) = @ARGV; $strat //= 'big'; setp($prime) if $prime;
my ($tn, $td) = split m{/}, $ts; $td //= 1; my $tau = mfrac($tn, $td);
my $m; my @yv; srand(12345 + $K);
if ($mode eq 'sym') { $m = $K; for my $L (1..$K) { my @v = (0) x $m; $v[$L-1] = 1; $yv[$L] = \@v; } }
else { ($m) = $mode =~ /^r(\d+)$/ or die "mode\n"; for my $L (1..$K) { $yv[$L] = [ map { 1 + int(rand($P - 2)) } 1..$m ]; } }
my $t0 = time;
if ($strat eq 'cons') {
  my $S1 = kp_setup($tau, $K, \@yv, 'big'); my $S2 = kp_setup($tau, $K, \@yv, 'small'); my $S3 = kp_setup($tau, $K, \@yv, 'mix');
  my ($n, $bad) = (0, 0);
  for my $k (keys %{$S1->{Y}}) { $n++; my ($v1, $v2, $v3) = ($S1->{Y}{$k}, $S2->{Y}{$k}, $S3->{Y}{$k});
    for my $i (0..$m-1) { if ($v1->[$i] != $v2->[$i] || $v1->[$i] != $v3->[$i]) { $bad++; print "INCONSISTENT Y($k)\n" if $bad < 10; last; } } }
  printf "ref_t1 cons K=%d tau=%s p=%d: %d multisets, %d inconsistent (strategies big/small/mix), %ds\n", $K, $ts, $P, $n, $bad, time - $t0; exit; }
my $S = kp_setup($tau, $K, \@yv, $strat); my $Pi = $S->{Pi};
my $Pj = leaf_series($tau, $Pi, $K, $K + 1);
# P^T/|Aut T| series, built by adding the smallest part
my %PS; my @Ts = sort { scalar(@$a) <=> scalar(@$b) } @{$S->{Ts}};
for my $T (@Ts) { my @t = @$T; my $s = pop @t; my $mult = scalar grep { $_ == $s } @$T;
  my $base = @t ? $PS{keyof(@t)} : [1, (0) x $K]; my $pr = ser_mul($base, $Pj->{$s}, $K); my $im = mi($mult);
  $PS{keyof(@$T)} = [ map { mm($_, $im) } @$pr ]; }
my @rel = map { [ (0) x $m ] } 0..$K;
for my $T (@Ts) { my $k = keyof(@$T); my $ps = $PS{$k}; my $y = $S->{Y}{$k};
  for my $o (Dof(@$T)..$K) { my $c = $ps->[$o]; next unless $c; for my $i (0..$m-1) { $rel[$o][$i] = ($rel[$o][$i] + $c*$y->[$i]) % $P if $y->[$i]; } } }
my $kap = kappa_list($tau, $K);
my ($n, $bad) = (0, 0); my @badk;
for my $k (1..$K) { my @pred = (0) x $m;
  for my $L (1..$k) { my $N = $L + 1; my $d = $k - $L;
    my $beta = mm($kap->[$k], mm(($d % 2) ? $P - 1 : 1, md(mm(md(1 + $tau), gbin(mm($N, $tau), $d)) - mm($tau, gbin(md(mm($N, md(1 + $tau)) + $d - 1), $d)))));
    for my $i (0..$m-1) { $pred[$i] = ($pred[$i] + $beta*$yv[$L][$i]) % $P; } }
  my $ok = 1; for my $i (0..$m-1) { $n++; if ($pred[$i] != $rel[$k][$i]) { $bad++; $ok = 0; } } push @badk, $k unless $ok; }
my $nz = 0; for my $k (1..$K) { $nz++ if grep { $_ } @{$rel[$k]}; }
printf "ref_t1 K=%d tau=%s mode=%s strat=%s p=%d: %d multisets, %d coefficients compared, %d mismatches%s; non-zero relation orders %d; %ds\n",
  $K, $ts, $mode, $strat, $P, scalar(@Ts), $n, $bad, (@badk ? " (orders @badk)" : ""), $nz, time - $t0;
