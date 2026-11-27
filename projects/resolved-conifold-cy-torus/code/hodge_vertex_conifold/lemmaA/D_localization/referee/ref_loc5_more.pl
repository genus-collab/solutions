# ref_loc5_more.pl -- larger cases of ref_loc5.pl (same code), torus a1 = e1, alpha0 = 1, and generic T_CY.
use strict; use warnings; use FindBin;
{ local @ARGV = (); no warnings; }
our $NORUN = 1;
do "$FindBin::Bin/ref_loc5_core.pl" or die $@ || $!;
my ($e1, $e2) = (Q(2,3), Q(5,4)); my $bad = 0;
for my $T (["a1 = e1, alpha0 = 1", [Q(1), $e1, $e2 - 1]], ["generic T_CY, alpha0 = 1, a1 = 1/3", [Q(1), Q(1,3), $e1 + $e2 - 1 - Q(1,3)]]) {
  my ($name, $w) = @$T; my @W = (@$w, $e1, $e2); print "torus: $name\n";
  for my $cs ([2,3], [3,2], [4,1]) { my ($g, $d) = @$cs; my $t0 = time; my $got = GW5($g, $d, \@W);
    my $pred = P_g($g, $e1, $e2) * Q($d)**(2*$g - 3) / ($e1 * $e2); my $ok = ($got == $pred); $bad++ unless $ok;
    printf "   (g,d) = (%d,%d): localization %s, prediction %s  %s [%d s]\n", $g, $d, $got, $pred, $ok ? "ok" : "MISMATCH", time - $t0; } }
print $bad ? "SOME MISMATCH\n" : "ALL OK\n";
