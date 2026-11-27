#!/usr/bin/perl
# ref_special.pl -- referee's check of Lemma 8.1 (special-framing lemma) with the referee's own two-leaf values:
# at tau0 = -k0/q (gcd 1, 1 <= k0 < q <= 2G-2), the row space of the two-leaf regularity conditions
#   { Y(m,n)(tau0) : q | n, m >= 1, m+n <= 2G-1 }   (Y from refkp's propagation; m = 1 gives n y_{n-1})
# equals the row space of { [y^N] (1-y)^(k0-1) (1-(1+tau0) y) F(y) : q <= N <= 2G-2 },  F(y) = Sum_c c y_{c-1} u(y)^c, u = y (1-y)^tau0.
# Usage: perl ref_special.pl Gmin Gmax   (mod p)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use refkp;
my ($Gmin, $Gmax) = @ARGV; set_mode($ENV{MODE} // 'mod');
sub rank_rows { my @M = map { [ @$_ ] } @{$_[0]}; my $nc = $_[1]; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { unless (IZ($M[$i][$c])) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = IV($M[$rk][$c]); $M[$rk] = [ map { ML($_, $inv) } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || IZ($M[$i][$c]); my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = SB($M[$i][$j], ML($f, $M[$rk][$j])); } } $rk++; } return $rk; }
my ($n, $bad) = (0, 0);
for my $G ($Gmin..$Gmax) { my $K = 2*$G - 3;
  for my $q (2..2*$G-2) { for my $k0 (1..$q-1) { next unless gcd($k0, $q) == 1; my $tau = NG(FR($k0, $q));
    my $S = kp_propagate($tau, $K);
    my @A; for (my $nn = $q; $nn <= 2*$G-2; $nn += $q) { for my $m (1..2*$G-1-$nn) { if ($m == 1) { push @A, vscale(NN($nn), vunit($K, $nn - 2)); } else { push @A, $S->{Y}{mkey($m, $nn)}; } } }
    # f_N = [y^N] F as linear forms (index c-2 <-> y_{c-1})
    my @f; for my $N (0..2*$G-2) { my $v = vzero($K); for my $c (2..$N) { next if $c - 1 > $K; my $e = $N - $c; vaxpy($v, ML(NN($c), ML(NN(($e % 2) ? -1 : 1), binom_gen(ML(NN($c), $tau), $e))), vunit($K, $c - 2)); } $f[$N] = $v; }
    # w(y) = (1-y)^(k0-1) (1 - (1+tau) y)
    my @w = (ONE()); for (1..$k0-1) { my @nw = map { ZERO() } 0..@w; for my $i (0..$#w) { $nw[$i] = AD($nw[$i], $w[$i]); $nw[$i+1] = SB($nw[$i+1], $w[$i]); } @w = @nw; }
    { my @nw = map { ZERO() } 0..@w; for my $i (0..$#w) { $nw[$i] = AD($nw[$i], $w[$i]); $nw[$i+1] = SB($nw[$i+1], ML(AD(ONE(), $tau), $w[$i])); } @w = @nw; }
    my @B; for my $N ($q..2*$G-2) { my $v = vzero($K); for my $i (0..$#w) { my $k = $N - $i; next if $k < 0; vaxpy($v, $w[$i], $f[$k]); } push @B, $v; }
    my $ra = @A ? rank_rows(\@A, $K) : 0; my $rb = @B ? rank_rows(\@B, $K) : 0; my $rab = rank_rows([ @A, @B ], $K);
    $n++; my $ok = ($ra == $rb && $rb == $rab); $bad++ unless $ok; print "G=$G q=$q k0=$k0: ranks $ra $rb $rab DIFFERENT\n" unless $ok; } }
  print "G=$G: framings so far $n, differing $bad\n"; }
print "ref_special G=$Gmin..$Gmax: $n framings, $bad differ\n";
