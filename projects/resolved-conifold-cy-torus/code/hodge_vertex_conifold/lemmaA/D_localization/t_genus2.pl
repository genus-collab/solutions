# t_genus2.pl -- check of xi_lib.pl against the proved genus-2 closed formula (proof_M_residue.md s.2)
#   5760 sigma^2 Xi_2(j) = (j-1)(j-2)[7 - (12j-28) sigma + (5j^2-26j+24) sigma^2],  sigma = tau(1+tau)
use strict; use warnings; use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }
my $tau = Q(2,7);
for my $R ([1],[2],[3],[0],[1,1],[2,2],[2,1]) { print "h=2 R=(@$R): Xi = ", xi(2, $tau, @$R), "\n"; }
for my $j (3,4,5) { my $sig = $tau*(1+$tau); my $pred = ($j-1)*($j-2)*(7 - (12*$j-28)*$sig + (5*$j*$j-26*$j+24)*$sig**2)/(5760*$sig**2);
  my $got = xi(2,$tau,$j); print "j=$j closed formula=$pred  engine=$got  ", ($pred == $got ? "ok" : "MISMATCH"), "\n"; }
