#!/usr/bin/perl
# ref_dump.pl -- dump Y(T) (linear forms in y_1..y_K) from the referee propagation, exact.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use refkp;
my ($K, $ts) = @ARGV; set_mode('rat'); my ($a, $b) = split m{/}, $ts; $b //= 1; my $S = kp_propagate(FR($a, $b), $K);
for my $k (sort keys %{$S->{Y}}) { print "Y $k : ", join(" ", map { "$_" } @{$S->{Y}{$k}}), "\n"; }
