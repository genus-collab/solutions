# slice2.pl h m : the slice Xi_h(m, j) of the two-leaf function (sigma-form sigma^h Xi), divided by its forced zeros
#   j = 0, 1, ..., 2h-1-m  (string, dilaton and Lemma A), and compared with the one-leaf bracket br_h(j) = sigma^h Xi_h(j)/prod_{l=1}^{2h-2}(j-l).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly;
sub R { Math::BigRat->new(shift) }
my ($h, $m) = @ARGV; $m //= 2;
my $C2 = CPoly::load($FindBin::Bin, $h, 2); my $S2 = CPoly::sigma_form($C2, $h);
my $C1 = CPoly::load($FindBin::Bin, $h, 1); my $S1 = CPoly::sigma_form($C1, $h);
my $sl = CPoly::subst_poly($S2, { j1 => Poly::const($m) });          # poly in j2, s
# divide by (j2 - l), l = 0..2h-1-m
sub divlin { my ($p, $v, $a) = @_;   # p / (v - a), exact; returns quotient
  my $d = Poly::degree_in($p, $v); my @c = map { Poly::coeff_of($p, $v, $_) } 0..$d;   # coefficients (polys in other vars)
  my @q; my $carry = {}; for (my $k = $d; $k >= 1; $k--) { $carry = Poly::add($c[$k], Poly::scale($carry, $a)); $q[$k-1] = $carry; }
  my $rem = Poly::add($c[0], Poly::scale($carry, $a)); die "not divisible by ($v - $a)" unless Poly::is_zero($rem);
  my $res = {}; for my $k (0..$#q) { $res = Poly::add($res, Poly::mul($q[$k], Poly::pow(Poly::var($v), $k))); } return $res; }
my $Q = $sl; for my $l (0..2*$h-1-$m) { $Q = divlin($Q, 'j2', $l); }
print "genus $h: sigma^h Xi_h($m, j) / prod_{l=0}^{2h-1-$m} (j-l) = Q(j), degree in j ", Poly::degree_in($Q, 'j2'), "\n";
my $br = $S1; for my $l (1..2*$h-2) { $br = divlin($br, 'j1', $l); }
# print Q and br in (j, s)
sub show { my ($p, $v) = @_; for my $k (0..Poly::degree_in($p, $v)) { my $c = Poly::coeff_of($p, $v, $k); my @t;
    for my $sp (0..$h) { my $x = Poly::const_val(Poly::coeff_of($c, 's', $sp)); push @t, "$x*s^$sp" if $x != 0; } print "   $v^$k: ", join(' + ', @t), "\n"; } }
print " Q(j) coefficients (j = j2):\n"; show($Q, 'j2');
print " br(j) coefficients (j = j1):\n"; show($br, 'j1');
# test: Q(j) - c * br(j+1) for the best constant c read off the top coefficient? print Q(j)/br(j+1) at a few j, s
for my $jv (2, 5, 9) { for my $sv (R(1), R(2)) { my $q = Poly::const_val(Poly::subst($Q, { j2 => R($jv), s => $sv })); my $b1 = Poly::const_val(Poly::subst($br, { j1 => R($jv + 1), s => $sv })); my $b0 = Poly::const_val(Poly::subst($br, { j1 => R($jv), s => $sv }));
  print "  j = $jv, sigma = $sv: Q = $q, br(j+1) = $b1, br(j) = $b0, Q/br(j+1) = ", ($b1 != 0 ? $q/$b1 : 'inf'), "\n"; } }
