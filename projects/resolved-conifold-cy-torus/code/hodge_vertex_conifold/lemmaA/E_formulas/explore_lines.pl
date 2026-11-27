# explore_lines.pl: for r = 2, which hyperplanes J = l (l integer) divide each component Xi^(k)_h(j1, j2)?
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly;
for my $h (2..6) { my $C = CPoly::load($FindBin::Bin, $h, 2);
  for my $k (0..$h) { my @ls;
    for my $l (-3..(3*$h+2)) { my $q = CPoly::subst_poly($C->[$k], { j2 => Poly::sub_(Poly::const($l), Poly::var('j1')) }); push @ls, $l if Poly::is_zero($q); }
    print "genus $h, k = $k: Xi^(k)(j1, l - j1) == 0 for l in {@ls}\n"; } }
