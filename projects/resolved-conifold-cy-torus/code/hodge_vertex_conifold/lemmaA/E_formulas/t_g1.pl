use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt;
Hodge::setup(1);
for my $j (0..4) { my $v = HodgePt::WP([], [$j]); print "genus 1, j = $j: components ", join(", ", @$v), "   expected Xi_1 = (j - 2 - 1/sigma)/24 -> k0 = ", ($j-2)/24, ", k1 = 1/24\n"; }
for my $J ([2,2],[2,3],[1,3],[0,3]) { my $v = HodgePt::WP([], $J); print "genus 1, R = (@$J): ", join(", ", @$v), "\n"; }
