#!/usr/bin/perl
# beta_kappa.pl -- divide f_{k,L}(sigma) by kappa_k = B_k/(k k!) and print the sigma-coefficients, grouped by depth d = k - L.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my @Bn = (R(1)); for my $m (1..30) { my $s = R(0); for my $k (0..$m-1) { $s = $s + binom($m+1, $k) * $Bn[$k]; } $Bn[$m] = -$s/($m+1); }
my %F; for (`perl "$FindBin::Bin/beta_sigma.pl"`) { next unless /^f_\{(\d+),(\d+)\} = (theta\*)?\[ (.*) \]/; $F{"$1,$2"} = [ map { R($_) } split /, /, $4 ]; }
for my $d (2..12) { print "depth $d:\n";
  for my $k (sort { $a <=> $b } map { (split /,/)[0] } grep { my ($k,$L) = split /,/; $k - $L == $d } keys %F) { my $L = $k - $d;
    my $kap = $Bn[$k] / ($k * fact($k)); my @c = map { $_ / $kap } @{$F{"$k,$L"}};
    printf "  k=%2d L=%2d : %s\n", $k, $L, join(', ', map { "$_" } @c); } }
