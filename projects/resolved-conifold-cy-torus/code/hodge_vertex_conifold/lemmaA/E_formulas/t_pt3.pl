use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt;
Hodge::setup(2);
sub sumW { my ($beta, $j) = @_; my $s = [Math::BigRat->new(0)]; for my $a (0..6) { my $v = Hodge::W($beta, [$a]); $s = Hodge::vadd($s, Hodge::vscale($v, Math::BigRat->new($j) ** $a)); } return $s; }
for my $beta ([], [1], [2], [1,1]) { for my $j (0, 1, 2, 3) {
  my $a = HodgePt::WP($beta, [$j]); my $b = sumW($beta, $j);
  print "beta=(@$beta) j=$j: WP ", join(",", @$a), "   W ", join(",", @$b), "\n"; } }
