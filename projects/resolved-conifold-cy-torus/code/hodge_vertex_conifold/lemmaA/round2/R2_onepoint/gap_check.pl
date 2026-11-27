#!/usr/bin/perl
# gap_check.pl ENGINE_FILE h : the r = 0 gap formula implied by (6.1) at j = 0 and the string/dilaton relation Xi_h(0) = (2-2h) F_h,
#   F_h = Int_{Mbar_h} K Omega = -(2h-3)! [lambda^{2h}] S(e1 lambda) S(e2 lambda) / sigma^h,  e1 = -tau, e2 = 1 + tau,
# compared with Xi_h(0) read from Route E's engine table engine_h<h>_r1.txt (j^0 coefficient), at four values of tau.
use strict; use warnings; use Math::BigRat; sub Q { Math::BigRat->new(@_) }
my ($file, $h) = @ARGV; my %C; open my $fh, '<', $file or die "no $file"; while (my $l = <$fh>) { chomp $l; my @f = split /\s*\|\s*/, $l; next unless @f >= 4; for my $t (split /\s*;\s*/, $f[3]) { next unless $t =~ /^(\d+):(\S+)$/; $C{$f[2]}{$1} = Q($2); } }
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; $r }
sub bern { my $n = shift; my @B = (Q(1)); for my $m (1..$n) { my $s = Q(0); for my $k (0..$m-1) { $s += fact($m+1)/(fact($k)*fact($m+1-$k))*$B[$k]; } $B[$m] = -$s/($m+1); } $B[$n] }
for my $tau (map { Q($_) } ("1", "2", "-3/7", "5/2")) { my $sig = $tau*(1+$tau); my ($e1,$e2) = (-$tau, 1+$tau);
  my $xi0 = Q(0); for my $k (keys %C) { $xi0 += ($C{$k}{0} // Q(0)) * (-Q(1)/$sig)**$k; }
  my @lc = map { Q(0) } 0..2*$h; for my $n (1..$h) { $lc[2*$n] = -bern(2*$n)/(2*$n*fact(2*$n)) * ($e1**(2*$n) + $e2**(2*$n)); }
  my @ex = (Q(1)); for my $n (1..2*$h) { my $a = Q(0); $a += $_*$lc[$_]*$ex[$n-$_] for 1..$n; $ex[$n] = $a/$n; }
  my $gap = -fact(2*$h-3) * $ex[2*$h] / $sig**$h;
  printf "h=%d tau=%s: Xi_h(0)/(2-2h) = %s ; gap formula = %s ; %s\n", $h, $tau, $xi0/(2-2*$h), $gap, ($xi0/(2-2*$h) == $gap ? "OK" : "MISMATCH"); }
