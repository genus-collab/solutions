# hodge_known.pl : checks of the Hodge engine against closed formulas (genus 2..6):
#  Faber: Int_{Mbar_g} lambda_g lambda_{g-1} lambda_{g-2} = |B_{2g-2}| |B_{2g}| / (2 (2g-2)! (2g-2) (2g))
#  Getzler-Pandharipande: Int_{Mbar_{g,1}} psi^{g-1} lambda_g lambda_{g-1} = |B_{2g}| / (2^{2g-1} (2g-1)!! (2g))
#  lambda_g formula: Int_{Mbar_{g,1}} psi^{2g-2} lambda_g = b_g = (2^{2g-1}-1)|B_{2g}| / (2^{2g-1} (2g)!)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge;
sub R { Math::BigRat->new(shift) }
my %B = (2 => R("1/6"), 4 => R("1/30"), 6 => R("1/42"), 8 => R("1/30"), 10 => R("5/66"), 12 => R("691/2730"));
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x } sub df { my $n = shift; my $x = R(1); for (my $k = $n; $k > 1; $k -= 2) { $x *= $k } $x }
# integrate a lambda-monomial (list of indices) times psi^a, via the ch-expansion of the lambdas
sub lamint { my ($g, $lams, $a) = @_; my @lam = Hodge::lambdas($g, 3*$g - 3 + scalar(@$a)); my $p = { '0,0,0' => [R(1)] };
  for my $k (@$lams) { $p = Hodge::chmul($p, $lam[$k], 3*$g - 3 + scalar(@$a)); }
  my $s = R(0); for my $mk (keys %$p) { my @m = split /,/, $mk; my $c = $p->{$mk}[0]; next if $c == 0; $s += $c * Hodge::HK($g, \@m, [], $a); } return $s; }
for my $g (2..6) {
  my $f1 = lamint($g, [$g, $g-1, $g-2], []); my $e1 = $B{2*$g-2} * $B{2*$g} / (2 * fact(2*$g-2) * (2*$g-2) * (2*$g));
  my $f2 = lamint($g, [$g, $g-1], [$g-1]); my $e2 = $B{2*$g} / (R(2)**(2*$g-1) * df(2*$g-1) * (2*$g));
  my $f3 = lamint($g, [$g], [2*$g-2]); my $e3 = (R(2)**(2*$g-1) - 1) * $B{2*$g} / (R(2)**(2*$g-1) * fact(2*$g));
  printf "genus %d: l_g l_{g-1} l_{g-2} %s | psi^{g-1} l_g l_{g-1} %s | psi^{2g-2} l_g %s\n", $g, ($f1 == $e1 ? "OK" : "FAIL $f1 vs $e1"), ($f2 == $e2 ? "OK" : "FAIL $f2 vs $e2"), ($f3 == $e3 ? "OK" : "FAIL $f3 vs $e3"); }
