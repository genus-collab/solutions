# quick sanity test of RefTR: tau = 2, compare with known values (omega_{0,3} = -sigma at the double-pole basis,
# Xi_1(j) = (j-2)/24 - 1/(24 sigma), Xi_2(3) = -73/20736 at tau = 2).
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use RefMod; use RefTR;
my @P = prime_list(1); setp($P[0]);
my $C = RefTR::new_curve(2,1,6*2+2*1-4); RefTR::run($C,2,1,verbose=>1);
print "omega_03: ", join(" ", map { "$_=>".$C->{W}{"0,3"}{$_} } keys %{$C->{W}{"0,3"}}), "  (-6 mod p = ".negm(6).")\n";
for my $j (1..3) { my $v = RefTR::xi_p2($C,1,$j); my $e = subm(ratm($j-2,24), ratm(1,144)); print "Xi_1($j): ", ($v == $e ? "ok" : "BAD $v vs $e"), "\n"; }
my $v = RefTR::xi_p2($C,2,3); print "Xi_2(3): ", ($v == ratm(-73,20736) ? "ok" : "BAD"), "\n";
for my $j (1..2) { print "Xi_2($j) = ", RefTR::xi_p2($C,2,$j), "\n"; }
