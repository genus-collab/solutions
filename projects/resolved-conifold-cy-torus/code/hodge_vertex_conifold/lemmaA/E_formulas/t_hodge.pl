use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Hodge; use Time::HiRes qw(time);
my $t0 = time;
print "lambda1 on M11: ", Hodge::HK(1,[1,0,0],[],[0]), " (1/24)\n";
print "lambda1^3 on M2: ", Hodge::HK(2,[3,0,0],[],[]), " (1/2880)\n";
print "ch3 on M2: ", Hodge::HK(2,[0,1,0],[],[]), " (-1/34560)\n";
print "lambda1 psi^3 on M21: ", Hodge::HK(2,[1,0,0],[],[3]), " (1/480)\n";
print "lambda1^6 on M3: ", Hodge::HK(3,[6,0,0],[],[]), " (1/90720)  t=", time-$t0, "\n";
print "lambda1^9 on M4: ", Hodge::HK(4,[9,0,0],[],[]), " (1/113400)  t=", time-$t0, "\n";
print "memo sizes ", scalar(keys %Hodge::memo), " ", scalar(keys %Hodge::memo0), " ", scalar(keys %WK::memo), "\n";
