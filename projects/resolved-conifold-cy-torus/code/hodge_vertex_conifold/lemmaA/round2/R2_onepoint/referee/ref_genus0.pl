#!/usr/bin/perl
# genus-0 check of formal Theorem 1 + P2: Xi_0(R) = Int_{Mbar_{0,r}} K prod 1/(1 - j_i psi_i) = (J-1)(J-2)...(J-r+3), mod p
use strict; use warnings; my $P = 2147483647; my ($tot, $bad) = (0, 0);
sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
for my $t ("1/2", "-2/7", "3", "5/3") { for my $cfg ([3, 6], [4, 5], [5, 4], [6, 3]) { my ($r, $jm) = @$cfg;
  for (`perl ref_xy_mod.pl $t 0 $r $jm -`) { next unless /mu=\(([\d,]+)\): ref (\d+)/; my @m = split /,/, $1; my $v = $2; my $J = 0; $J += $_ for @m;
    my $p = 1; $p = ($p * ($J - $_)) % $P for 1..($r - 3); $tot++; if ($p != $v) { $bad++; print "MISMATCH tau=$t mu=(@m)\n"; } } } }
print "genus 0, r = 3..6, 4 framings: compared $tot with (J-1)...(J-r+3), mismatches $bad\n";
