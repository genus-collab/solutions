#!/usr/bin/perl
# ref_forcing.pl -- referee 2: direct test of the special-framing lemma (s.8) and of the forcing theorem (s.9) at genus G,
# WITHOUT using T1 or Lemma 8.1: the polar relations are computed from the full KP propagation (rk.pm) on the leaf curve,
# the two-leaf values from the same propagation (Theorem D, R = empty).
# At tau0 = -k0/q (gcd 1, 2k0 < q), level D odd, unknowns y_D..y_{2G-3} (y_L = 0 for L < D):
#   (TL) two-leaf/one-leaf regularity: Y(m,n)(tau0) = 0 if q | m or q | n (m + n <= 2G-1), y_{c-1}(tau0) = 0 if q | c (c <= 2G-2);
#   (PO) even polar relations: relation_k(tau0) = 0, k even, 2 <= k <= 2G-4 (Route D input (I1));
#   question: is y_D forced to 0?
# Also compares, at each framing, the row space of (TL) with the author's p-conditions [y^N] p = 0, q <= N <= 2G-2,
#   p = (1-y)^(k0-1)(1-(1+tau0)y) F(y), F = sum_c c y_{c-1} u(y)^c, u = y(1-y)^tau0   (Lemma 8.1), with D = 1 (all unknowns).
# Usage: perl ref_forcing.pl G [all]     ('all': also framings outside the hypothesis range, reported separately)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rk;
my ($G, $all) = @ARGV; my $Lmax = 2*$G - 3;
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
sub rank { my ($rows, $nc) = @_; my @M = map { [ @$_ ] } @$rows; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { if ($M[$i][$c]) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = mi($M[$rk][$c]); $M[$rk] = [ map { mm($_, $inv) } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || !$M[$i][$c]; my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = md($M[$i][$j] - mm($f, $M[$rk][$j])); } } $rk++; } return $rk; }
sub polar_rows { my ($S, $tau, $K, $m) = @_; my $Pj = leaf_series($tau, $S->{Pi}, $K, $K + 1); my %PS; my @Ts = sort { scalar(@$a) <=> scalar(@$b) } @{$S->{Ts}};
  for my $T (@Ts) { my @t = @$T; my $s = pop @t; my $mult = scalar grep { $_ == $s } @$T; my $base = @t ? $PS{keyof(@t)} : [1, (0) x $K];
    my $pr = ser_mul($base, $Pj->{$s}, $K); my $im = mi($mult); $PS{keyof(@$T)} = [ map { mm($_, $im) } @$pr ]; }
  my @rel = map { [ (0) x $m ] } 0..$K;
  for my $T (@Ts) { my $k = keyof(@$T); my $ps = $PS{$k}; my $y = $S->{Y}{$k};
    for my $o (Dof(@$T)..$K) { my $c = $ps->[$o]; next unless $c; for my $i (0..$m-1) { $rel[$o][$i] = ($rel[$o][$i] + $c*$y->[$i]) % $P if $y->[$i]; } } }
  return \@rel; }
my %stat; my $t0 = time; my @lemma; my ($lemtests, $lembad) = (0, 0);
for my $q (3..2*$G-2) { for my $k0 (1..int(($q-1)/2)) { next unless gcd($k0, $q) == 1; my $tau = mfrac(-$k0, $q);
  # ---- Lemma 8.1 comparison with all unknowns y_1..y_Lmax
  { my $m = $Lmax; my @yv; for my $L (1..$Lmax) { my @v = (0) x $m; $v[$L-1] = 1; $yv[$L] = \@v; }
    my $S = kp_setup($tau, $Lmax, \@yv, 'big'); my @TL;
    for my $c (2..2*$G-2) { push @TL, $yv[$c-1] if $c % $q == 0; }
    for my $s (4..2*$G-1) { for my $n2 (2..int($s/2)) { my $m1 = $s - $n2; next unless $m1 % $q == 0 || $n2 % $q == 0; push @TL, $S->{Y}{keyof($m1, $n2)}; } }
    # p-conditions
    my @w = (1); for (1..$k0-1) { my @nw = (0) x (@w + 1); for my $i (0..$#w) { $nw[$i] = md($nw[$i] + $w[$i]); $nw[$i+1] = md($nw[$i+1] - $w[$i]); } @w = @nw; }
    { my @nw = (0) x (@w + 1); for my $i (0..$#w) { $nw[$i] = md($nw[$i] + $w[$i]); $nw[$i+1] = md($nw[$i+1] - mm(md(1 + $tau), $w[$i])); } @w = @nw; }
    my @f; for my $N (1..2*$G-2) { my @r = (0) x $m; for my $c (2..$N) { next if $c - 1 > $Lmax; $r[$c-2] = mm($c, mm(($N - $c) % 2 ? $P - 1 : 1, gbin(mm($c, $tau), $N - $c))); } $f[$N] = \@r; }
    my @PC; for my $N ($q..2*$G-2) { my @r = (0) x $m; for my $i (0..$#w) { my $k = $N - $i; next if $k < 1; for my $j (0..$m-1) { $r[$j] = md($r[$j] + mm($w[$i], $f[$k][$j])); } } push @PC, \@r; }
    my ($r1, $r2, $r12) = (rank(\@TL, $m), rank(\@PC, $m), rank([@TL, @PC], $m)); $lemtests++;
    unless ($r1 == $r2 && $r2 == $r12) { $lembad++; print "LEMMA 8.1 DIFFERS: G=$G q=$q k0=$k0 ranks $r1 $r2 $r12\n"; } }
  # ---- forcing at each odd D
  for (my $D = 3; $D <= $Lmax; $D += 2) { next if $q < $D + 2 && !$all;
    my $hyp = ($k0 <= $D + 1 && 3*($q - $k0) < 2*$G && $q >= $D + 2) ? 'in' : 'out'; next if $hyp eq 'out' && !$all;
    next if ($D + 1) % $q == 0;   # Pi_{D+1}(tau0) = 0: y_D vanishes trivially, no information on Xi_G(D+1)
    my $m = $Lmax - $D + 1; my @yv; for my $L (1..$Lmax) { my @v = (0) x $m; $v[$L-$D] = 1 if $L >= $D; $yv[$L] = \@v; }
    my $S = kp_setup($tau, $Lmax, \@yv, 'big'); my @rows;
    for my $c (2..2*$G-2) { push @rows, $yv[$c-1] if $c % $q == 0; }
    for my $s (4..2*$G-1) { for my $n2 (2..int($s/2)) { my $m1 = $s - $n2; next unless $m1 % $q == 0 || $n2 % $q == 0; push @rows, $S->{Y}{keyof($m1, $n2)}; } }
    my $rel = polar_rows($S, $tau, $Lmax, $m); for (my $k = 2; $k <= 2*$G-4; $k += 2) { push @rows, $rel->[$k]; }
    # odd polar orders must vanish identically (kappa_k = 0) -- record if not
    for (my $k = 3; $k <= $Lmax; $k += 2) { if (grep { $_ } @{$rel->[$k]}) { $stat{oddnonzero}++; } }
    my $r1 = rank(\@rows, $m); my @rows2 = map { [ @{$_}[1..$m-1] ] } @rows; my $r2 = $m > 1 ? rank(\@rows2, $m - 1) : 0;
    my $forced = ($r1 - $r2 == 1) ? 1 : 0; $stat{"$hyp,tot"}++; $stat{"$hyp,forced"} += $forced;
    print "NOT FORCED (inside hypotheses): G=$G D=$D q=$q k0=$k0\n" if $hyp eq 'in' && !$forced; } } }
printf "ref_forcing G=%d p=%d: Lemma 8.1 row spaces at %d framings, %d differ; forcing inside hypotheses: %d of %d forced%s; odd polar orders non-zero: %d; %ds\n",
  $G, $P, $lemtests, $lembad, $stat{'in,forced'} // 0, $stat{'in,tot'} // 0,
  ($all ? sprintf("; outside: %d of %d forced", $stat{'out,forced'} // 0, $stat{'out,tot'} // 0) : ""), $stat{oddnonzero} // 0, time - $t0;
