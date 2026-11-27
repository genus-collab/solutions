# boundary2.pl h : first non-vanishing layer D = 2h-2 for two leaves: sigma^h Xi_h(1+a, 2h-1-a)/(2h-1)!, a = 0..h-1, as sigma-polys,
# compared with [lambda^{2h}] S1 S2 B^{2h}  (the value for a = 0 given by the string equation and (6.1)).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Comps; use SerS;
sub R { Math::BigRat->new(shift) }
my $h = shift; my $C = Comps::load("$FindBin::Bin/engine_h${h}_r2.txt");
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
my $ref = SerS::coeff(SerS::smul(SerS::S1S2(), SerS::spow(SerS::sexp(SerS::beta()), 2*$h)), $h);
print "genus $h: [l^2h] S1S2 B^{2h} = ", join(' ', map { "($ref->[$_])s^$_" } grep { $ref->[$_] != 0 } 0..$#$ref), "\n";
for my $a (0..$h-1) { my @p = (1 + $a, 2*$h - 1 - $a); my @s;
  for my $k (0..$h) { $s[$h - $k] = (($k % 2) ? -1 : 1) * Comps::comp_at($C, $h, 2, $k, \@p) / fact(2*$h - 1); }
  print "  a=$a (R = @p): ", join(' ', map { "($s[$_])s^$_" } grep { $s[$_] != 0 } 0..$#s), "\n";
  my @d = map { $s[$_] - ($ref->[$_] // 0) } 0..$h; print "      minus ref: ", join(' ', map { "($d[$_])s^$_" } grep { $d[$_] != 0 } 0..$#d), "\n"; }
