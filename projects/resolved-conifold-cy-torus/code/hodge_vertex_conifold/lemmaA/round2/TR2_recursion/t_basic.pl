# t_basic.pl: tests of TRtau.pm against known values.
# usage: perl t_basic.pl tau G [nmax]
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use TRtau qw(R new_curve run_tr xi_tr);
my $tau = shift // 1; my $G = shift // 3; my $NM = shift // 1;
my $t0=time;
my $C = new_curve($tau, 3*$G-3+$NM);
print "tau=$tau sigma=".$C->{sig}." setup ".(time-$t0)."s\n";
print "omega03: ", join(" ", map { "$_=>".$C->{W}{"0,3"}{$_} } keys %{$C->{W}{"0,3"}}), "\n";
print "omega11: ", join(" ", map { "$_=>".$C->{W}{"1,1"}{$_} } sort keys %{$C->{W}{"1,1"}}), "\n";
run_tr($C,$G,$NM,verbose=>1);
print "recursion ".(time-$t0)."s\n";
for my $g (2..$G){ my @out; for my $j (1..2*$g+1){ push @out, "Xi_$g($j)=".xi_tr($C,$g,$j); } print join("  ",@out),"\n"; }
