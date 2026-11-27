#!/usr/bin/perl
# a3_zcount.pl -- number Z(D) of zeros (with multiplicity, at sigma != 0) forced on sigma^G Xi_G(D+1) by the KP merge formula
#   Xi_G(T) = (J-1)...(J-r+1) * Pi_{D+1}/prod_{t in T} Pi_t * Xi_G(D+1)   (T at level D, all lower levels vanishing),
# i.e. the maximal pole order of Pi_{D+1}/prod Pi_t over T, summed over the special framings tau = -c/a (sigma = -c(a-c)/a^2).
# ord_{tau=-c/a} Pi_m = [a | m] (gcd(c,a) = 1);  at a = 2 (tau = -1/2) sigma-order = tau-order/2.
use strict; use warnings;
sub gcd { my ($a, $b) = @_; ($a, $b) = ($b, $a % $b) while $b; return $a; }
sub phi { my $n = shift; my $c = 0; for my $k (1..$n) { $c++ if gcd($k, $n) == 1; } return $c; }
print "  D   t=D+1   Z(D)   [contributions a:order x #sigma]\n";
for my $D (1..22) { my $Z = int($D/2); my @c = ("2:" . int($D/2) . "x1");
  for my $a (3..$D+1) { my $ord = int($D/($a-1)) - ((($D+1) % $a) == 0 ? 1 : 0); next if $ord <= 0; my $ns = phi($a)/2; $Z += $ord * $ns; push @c, "$a:${ord}x$ns"; }
  printf "%3d %6d %6d   %s\n", $D, $D+1, $Z, join(' ', @c); }
