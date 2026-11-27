#!/usr/bin/perl
# t_forcing.pl -- check of the FORCING THEOREM on its whole hypothesis range for G <= Gmax:
#   tau0 = -k0/q, gcd(k0,q) = 1, 1 <= k0 <= D+1, 2 k0 < q, q >= D+2, 3(q-k0) < 2G, D odd in [3, 2G-3]:
#   the odd-order closed-form relations (T1 at tau0, F = p/((1-y)^(k0-1)(1-(1+tau0)y)), deg p <= q-1, p = O(y^(D+1)))
#   force p_{D+1} = 0.  Matrix: M_{N,n} = (q-k0) binom(k0+N-n-1, N-n) + k0 binom(k0-1-n, N-n), odd N <= 2G-3, n = D+1..q-1.
# Also reports, for information, tuples OUTSIDE the hypotheses (k0 > D+1 or 3(q-k0) >= 2G) where forcing still holds.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($Gmax) = @ARGV;
sub gb { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
sub gcd { my ($x, $y) = @_; ($x, $y) = ($y, $x % $y) while $y; return $x; }
sub rank { my @M = map { [ @$_ ] } @{$_[0]}; my $nc = $_[1]; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { unless ($M[$i][$c]->is_zero) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = R(1)/$M[$rk][$c]; $M[$rk] = [ map { $_ * $inv } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || $M[$i][$c]->is_zero; my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = $M[$i][$j] - $f*$M[$rk][$j]; } } $rk++; } return $rk; }
my ($in, $inbad, $out, $outforced) = (0, 0, 0, 0);
for my $G (3..$Gmax) { for (my $D = 3; $D <= 2*$G-3; $D += 2) { for my $q ($D+2 .. 2*$G-2) { for my $k0 (1 .. int(($q-1)/2)) { next unless gcd($k0, $q) == 1;
    my $hyp = ($k0 <= $D+1 && 3*($q-$k0) < 2*$G) ? 1 : 0;
    next if !$hyp && !$ENV{OUTSIDE};
    my @cols = ($D+1 .. $q-1); my @rows;
    for (my $N = 1; $N <= 2*$G-3; $N += 2) { my @r = map { my $n = $_; $N < $n ? R(0) : R($q-$k0) * gb($k0 + $N - $n - 1, $N - $n) + R($k0) * gb($k0 - 1 - $n, $N - $n) } @cols; push @rows, \@r if grep { !$_->is_zero } @r; }
    my $nc = @cols; my $r1 = rank(\@rows, $nc); my @rows2 = map { [ @{$_}[1..$nc-1] ] } @rows; my $r2 = $nc > 1 ? rank(\@rows2, $nc-1) : 0;
    my $forced = ($r1 - $r2 == 1) ? 1 : 0;
    if ($hyp) { $in++; unless ($forced) { $inbad++; print "FAIL (inside hypotheses): G=$G D=$D q=$q k0=$k0\n"; } }
    else { $out++; $outforced += $forced; } } } } }
print "t_forcing G<=$Gmax: $in tuples inside the hypotheses, $inbad not forced", ($ENV{OUTSIDE} ? "; outside: $outforced of $out forced" : ""), "\n";
