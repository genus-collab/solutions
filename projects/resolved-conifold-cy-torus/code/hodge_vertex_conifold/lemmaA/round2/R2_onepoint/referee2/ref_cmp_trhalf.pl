#!/usr/bin/perl
# ref_cmp_trhalf.pl FILE : compare the exact EO one-point values on C_{-1/2} (output of status_trhalf.pl: [x^j]F_{g,1}/Pi_j = sigma^{1-g} Xi_g(j))
# with formula (6.1) evaluated by the referee's own implementation RefR2::xi61 (tau = -1/2, sigma = -1/4).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my $tau = Q("-1/2"); my $sig = $tau * (1 + $tau); my %cache; my ($N, $B) = (0, 0);
open my $f, '<', $ARGV[0] or die;
while (my $line = <$f>) { next unless $line =~ /g=(\d+) n=1: (.*)$/; my ($g, $rest) = ($1, $2); next if $g < 2; my ($n, $bad, $nz) = (0, 0, 0);
  while ($rest =~ /\((\d+)\)=(\S+)/g) { my ($j, $v) = ($1, Q($2)); $n++; $nz++ unless $v->is_zero;
    my $p = $sig**(1 - $g) * RefR2::xi61($g, $tau, $j, \%cache); if ($p != $v) { $bad++; print "  MISMATCH g=$g j=$j EO=$v (6.1)=$p\n"; } }
  $N += $n; $B += $bad; print "tau=-1/2 g=$g: $n odd j compared ($nz non-zero), mismatches $bad\n"; }
print "TOTAL $N, mismatches $B\n";
