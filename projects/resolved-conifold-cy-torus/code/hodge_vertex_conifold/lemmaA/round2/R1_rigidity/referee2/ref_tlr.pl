#!/usr/bin/perl
# ref_tlr.pl -- referee 2: spot check of the side claim "TLR nullity in the window G-1" (report s.12).
# Unknowns: c_{t,k}, sigma^G Xi_G(t) = sum_{k=1}^{W} c_{t,k} sigma^k, t = tmin..2G-2 (tmin = 3: Xi_G(2) = 0 by Theorem A; tmin = 2: Route A's set-up).
# Constraint: every two-leaf value sigma^G Xi_G(a,b) (a,b >= 2, a+b <= 2G-1), from Theorem D (R = empty) via rk.pm, lies in
# sigma Q[sigma]_{<= W-1}; imposed by Lagrange interpolation on W random nodes and evaluation at further random samples.
# Usage: perl ref_tlr.pl G W tmin [samples]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use rk;
my ($G, $W, $tmin, $S) = @ARGV; $S //= 2*$G + 8; srand(4242 + $G); my $L = 2*$G - 3;
my @ts = ($tmin..2*$G-2); my $nc = scalar(@ts) * $W; my %col; my $i = 0; for my $t (@ts) { for my $k (1..$W) { $col{"$t,$k"} = $i++; } }
sub rank { my ($rows, $nc) = @_; my @M = map { [ @$_ ] } @$rows; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { if ($M[$i][$c]) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = mi($M[$rk][$c]); $M[$rk] = [ map { mm($_, $inv) } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || !$M[$i][$c]; my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = md($M[$i][$j] - mm($f, $M[$rk][$j])); } } $rk++; } return $rk; }
my @samp;
for my $s (1..$S) { my $tau = 3 + int(rand($P - 9)); my $sig = mm($tau, md($tau + 1)); my $isg = mi($sig); my $sgG = mpow($isg, $G);
  my @Pi = (1, 1); for my $j (2..2*$G+2) { my $v = 1; for my $b (1..$j-1) { $v = mm($v, md($j*$tau + $b)); } $Pi[$j] = mm($v, ifct($j-1)); }
  my @yv; for my $Lv (1..$L) { my @v = (0) x $nc; my $t = $Lv + 1; if ($t >= $tmin) { my $pw = 1; for my $k (1..$W) { $pw = mm($pw, $sig); $v[$col{"$t,$k"}] = mm($Pi[$t], mm($sgG, $pw)); } } $yv[$Lv] = \@v; }
  my $St = kp_setup($tau, $L, \@yv, 'big'); my %two;
  for my $sum (4..2*$G-1) { for my $b (2..int($sum/2)) { my $a = $sum - $b; my $f = mm(mpow($sig, $G), mi(mm($Pi[$a], $Pi[$b]))); $two{"$a,$b"} = [ map { mm($f, $_) } @{$St->{Y}{keyof($a, $b)}} ]; } }
  push @samp, [$sig, \%two]; }
my @nodes = map { $samp[$_][0] } 0..$W-1; my @rows;
for my $t ($W..$S-1) { my $x = $samp[$t][0]; my @l;
  for my $s (0..$W-1) { my $v = mm($x, mi($nodes[$s])); for my $u (0..$W-1) { next if $u == $s; $v = mm($v, mm(md($x - $nodes[$u]), mi(md($nodes[$s] - $nodes[$u])))); } push @l, $v; }
  for my $k (keys %{$samp[$t][1]}) { my @r = @{$samp[$t][1]{$k}}; for my $s (0..$W-1) { my $v = $samp[$s][1]{$k}; for my $j (0..$nc-1) { $r[$j] = md($r[$j] - mm($l[$s], $v->[$j])); } } push @rows, \@r; } }
my $rk = rank(\@rows, $nc);
printf "ref_tlr G=%d window=%d t=%d..%d: %d unknowns, %d equations, rank %d, NULLITY %d\n", $G, $W, $tmin, 2*$G-2, $nc, scalar(@rows), $rk, $nc - $rk;
