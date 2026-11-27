use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod;
my $K = shift; my $tau = mfrac(2, 7); my $t0 = time; my $S = kp_setup($tau, 1, $K); my $t1 = time; my $rows = kp_rows($S); my $t2 = time;
print "K=$K: setup ", $t1-$t0, " s, rows ", $t2-$t1, " s, #T = ", scalar(@{$S->{allT}}), "\n";
