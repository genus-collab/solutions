# t_engine.pl: sanity tests of RefTR.pm (referee 2).
#  (a) omega_{0,3}, omega_{1,1} at tau = 2 (author: omega_{1,1} = (163/1728) e_0 - (1/4) e_1, e_k = (2k+1) xi_k);
#  (b) Xi_g(mu) at tau = 2 against the author's exact values; (c) genus 1, 2 one-leaf values against the Hodge formulas
#      sigma Xi_1(j) = ((j-2) sigma - 1)/24,  Xi_2(j) = (j-1)(j-2)[(5j^2-26j+24) + (28-12j)/sigma + 7/sigma^2]/5760;
#  (d) tau = 1 genus-5 values of the earlier referee: sigma^(1-g) Xi_5(9) = -2205/1048576.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefTR;
my $p = 2147483647; sub R { RefTR::rat($_[0],$_[1],$p) }
my $S = RefTR->new(2,1,3*4-2,$p); $S->run(4,1);
my $W11 = $S->{W}{"1,1"}; print "tau=2 omega_11 xi-coefs: k0=0 ", ($W11->{0}//0)==R(163,1728) ? "ok" : "BAD", "  k0=1 ", ($W11->{1}//0)==R(-3,4) ? "ok" : "BAD", "\n";
print "tau=2 omega_03 = ", $S->{W}{"0,3"}{"0,0,0"} == R(-6,1) ? "-6 ok" : "BAD", "\n";
my %want = ("2,3"=>[-73,20736], "3,5"=>[-245,82944], "4,7"=>[-182345,23887872]);
for my $k (sort keys %want) { my ($g,$mu)=split /,/,$k; my $v = $S->xi($g,$mu); print "tau=2 Xi_$g($mu): ", $v == R(@{$want{$k}}) ? "ok" : "BAD ($v)", "\n"; }
for my $g (2..4) { my @z; for my $mu (1..2*$g-2) { push @z, $S->xcoef($g,$mu) ? "nz" : "0"; } print "tau=2 g=$g [x^mu]F, mu=1..",2*$g-2,": @z\n"; }
# (c) genus 1,2 Hodge formulas at several tau
for my $tt ([2,1],[3,7],[-5,3],[5,2]) { my $T = RefTR->new($tt->[0],$tt->[1],5,$p); $T->run(2,1);
  my $sig = $T->{sig}; my $is = RefTR::inv($sig,$p); my $bad = 0;
  for my $j (1..7) { my $x1 = $T->xi(1,$j); next unless defined $x1; my $w1 = RefTR::md((($j-2)*$sig - 1) % $p * RefTR::inv(24,$p) % $p * $is, $p); $bad++ unless $x1 == $w1;
    my $x2 = $T->xi(2,$j); next unless defined $x2; my $br = RefTR::md(5*$j*$j - 26*$j + 24 + RefTR::md(28-12*$j,$p)*$is + 7*$is % $p * $is, $p);
    my $w2 = RefTR::md(($j-1)*($j-2) % $p * $br % $p * RefTR::inv(5760,$p), $p); $bad++ unless $x2 == $w2; }
  print "tau=$tt->[0]/$tt->[1]: genus 1,2 one-leaf vs Hodge formulas, j=1..7: ", $bad ? "$bad MISMATCHES" : "all agree", "\n"; }
# (d) tau = 1, genus 5
my $U = RefTR->new(1,1,3*5-2,$p); $U->run(5,1);
my $v = ($U->xcoef(5,9) * RefTR::inv(RefTR::gbin(17,8,$p),$p)) % $p; print "tau=1 [x^9]F_5/C(17,8) ", $v == R(-2205,1048576) ? "= -2205/1048576 ok" : "BAD", "\n";
