use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt;
Hodge::setup(2);
sub sumW { my ($beta, $j) = @_; my $s = [Math::BigRat->new(0)]; for my $a (0..6) { my $v = Hodge::W($beta, [$a]); $s = Hodge::vadd($s, Hodge::vscale($v, Math::BigRat->new($j) ** $a)); } return $s; }
for my $beta ([3], [1,2], [4], [1,3]) { my $a = HodgePt::WP($beta, [2]); my $b = sumW($beta, 2); print "beta=(@$beta) j=2: WP ", join(",", @$a), "   W ", join(",", @$b), "\n"; }
my $s = [Math::BigRat->new(0)]; for my $a1 (0..6) { for my $a2 (0..6) { $s = Hodge::vadd($s, Hodge::vscale(Hodge::W([], [$a1, $a2]), Math::BigRat->new(2) ** ($a1 + $a2))); } }
print "Sum_a 2^|a| W(a1,a2) = ", join(",", @$s), "\n";
print "W(3,0) = ", join(",", @{Hodge::W([], [3,0])}), "  W(0,3) = ", join(",", @{Hodge::W([], [0,3])}), "\n";
print "T(3) = ", join(",", @{Hodge::T(3)}), "  T(1,2) = ", join(",", @{Hodge::T(1,2)}), "\n";
