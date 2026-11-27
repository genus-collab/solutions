#!/usr/bin/perl
# beta_table.pl -- coefficients beta_{k,L}(tau) of the closed-form relations after KP substitution:
#   relation_k = sum_{L<=k} beta_{k,L} y_L,  y_L = Pi_{L+1} Xi_G(L+1)  (universal, independent of G).
# Implementation: run the MV-normalised KP propagation of kp_point.pl with D = 1 and print the relation rows.
use strict; use warnings; use FindBin;
my ($K, $tau) = @ARGV;
my $out = `perl "$FindBin::Bin/kp_point.pl" 1 $K $tau none 1 2>&1`;
print $out;
