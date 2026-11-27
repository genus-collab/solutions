use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use Comps;
Hodge::setup(2); my $C = Comps::load("$FindBin::Bin/engine_h2_r1.txt"); my $C2 = Comps::load("$FindBin::Bin/engine_h2_r2.txt");
for my $j (0..5) { my $v = HodgePt::WP([], [$j]); print "h=2 j=$j: WP ", join(",", @$v), "  poly ", join(",", map { Comps::comp_at($C, 2, 1, $_, [$j]) } 0..2), "\n"; }
for my $J ([0,0],[1,0],[2,0],[1,1],[2,1],[2,2],[3,0]) { my $v = HodgePt::WP([], $J); print "h=2 J=(@$J): WP ", join(",", @$v), "  poly ", join(",", map { Comps::comp_at($C2, 2, 2, $_, $J) } 0..2), "\n"; }
