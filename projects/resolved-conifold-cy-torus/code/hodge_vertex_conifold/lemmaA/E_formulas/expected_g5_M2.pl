# expected_g5_M2.pl : (main session, round-1 check M2) engine predictions (c-normalisation, as Mathematica expressions in x) for the genus-5 cross-check points
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt;
my $h = 5; Hodge::setup($h);
my @pts = ([4], [2,3], [3,4], [3,3,3], [4,6]);   # M2 points (predicted 0) and the control (4,6) at D = 2h-2
my @out;
for my $p (@pts) { my $v = HodgePt::WP([], $p); my $r = scalar @$p;
  # Xi^c(x) = (1+x)^{-(3h-3+2r)} Sum_k v^k Xi^(k),  v = -1/sigma = (1+x)^2/x
  my @t = map { "(" . ($v->[$_] // 0) . ")*((1+x)^2/x)^$_" } grep { ($v->[$_] // 0) != 0 } 0..$#$v;
  my $e = @t ? "(" . join(" + ", @t) . ")/(1+x)^" . (3*$h - 3 + 2*$r) : "0";
  push @out, "{" . join(",", @$p) . "} -> $e"; }
print "expected5 = {", join(",\n  ", @out), "};\n";
