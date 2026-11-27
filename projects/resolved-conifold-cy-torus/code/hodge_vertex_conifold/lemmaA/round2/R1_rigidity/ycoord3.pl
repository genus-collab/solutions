#!/usr/bin/perl
# ycoord3.pl -- is the KP solution "framing-0-like" in the y-coordinate?  For the solution with one-leaf data F(y) = y^n
# (F = V~(u(y)), u = y(1-y)^tau), compare the three-leaf function in y-coordinates,
#   W3~(y1,y2,y3) = sum_{t} Y(t1,t2,t3) prod u(y_i)^{t_i}   (ordered parts >= 1, unit parts by string),
# with the framing-0 merge value (n+1) y1 y2 y3 h_{n-1}(y1,y2,y3).  Prints the difference coefficients.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($n, $K, $taus) = @ARGV; $K //= $n + 4; $taus //= '2/7'; my $tau = R($taus);
# one-leaf data: c y_{c-1} = [u^c] y(u)^n = (n/c) [y^{c-n}] (1-y)^{-c tau}  -> y_{c-1} = (n/c^2) binom(c tau + c - n - 1, c - n)
sub gb { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
my @y; for my $c (2..$K+1) { $y[$c-1] = ($c < $n) ? R(0) : R($n) / R($c*$c) * gb($tau*$c + $c - $n - 1, $c - $n); }
# multi-leaf values via kp_point.pl machinery (exact): run it and capture Y forms -> re-implement minimal propagation here
my $out = `perl "$FindBin::Bin/kp_dump.pl" 1 $K $taus`;
my %Y; for (split /\n/, $out) { next unless /^Y (\S+) : (.*)$/; my $key = $1; my @f = map { R($_) } split ' ', $2; my $v = R(0); for my $i (0..$#f) { $v = $v + $f[$i] * $y[$i+1] if defined $y[$i+1]; } $Y{$key} = $v; }
sub Yval { my @T = sort { $b <=> $a } @_; my @Tp = grep { $_ >= 2 } @T; my $nb = @T - @Tp; return R(0) unless @Tp; my $J = 0; $J += $_ for @Tp; my $lv = $J - @Tp;
  return R(0) if $lv > $K; my $v = $Y{join(',', @Tp)}; return R(0) unless defined $v; my $f = R(1); $f = $f * ($J + $_) for 0..$nb-1; return $v * $f; }
# W3~ coefficients up to total y-degree M
my $M = $n + 2 + 2;  # compare degrees n+2 .. n+4
my %W;   # key "i,j,k" -> coefficient of y1^i y2^j y3^k
my %ucoef; sub ucoef { my ($t, $m) = @_; return $ucoef{"$t,$m"} //= do { R((-1)**$m) * gb($tau*$t, $m) }; }   # [y^{t+m}] u^t = [y^m](1-y)^{t tau}
for my $t1 (1..$M) { for my $t2 (1..$M) { for my $t3 (1..$M) { next if $t1 + $t2 + $t3 > $M; my $Yv = Yval($t1, $t2, $t3); next if $Yv->is_zero;
  for my $m1 (0..$M) { for my $m2 (0..$M) { for my $m3 (0..$M) { my ($i, $j, $k) = ($t1+$m1, $t2+$m2, $t3+$m3); next if $i + $j + $k > $M;
    $W{"$i,$j,$k"} = ($W{"$i,$j,$k"} // R(0)) + $Yv * ucoef($t1,$m1) * ucoef($t2,$m2) * ucoef($t3,$m3); } } } } } }
# framing-0: (n+1) y1y2y3 h_{n-1}: coefficient 1*(n+1) for every (i,j,k) >= 1 with i+j+k = n+2
for my $deg ($n+2..$M) { my @d; for my $i (1..$deg) { for my $j (1..$deg-$i) { my $k = $deg - $i - $j; next if $k < 1; next if $i < $j || $j < $k;
      my $w = $W{"$i,$j,$k"} // R(0); my $w0 = ($deg == $n + 2) ? R($n+1) : R(0); push @d, "($i,$j,$k): " . ($w - $w0); } }
  print "degree $deg: W3~ - framing0 = ", join('  ', @d), "\n"; }
