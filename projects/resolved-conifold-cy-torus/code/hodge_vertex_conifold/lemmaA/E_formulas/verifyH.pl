# verifyH.pl hmax : end-to-end check of the symmetric two-leaf structure
#   sigma^h Xi_h(j1,j2) = prod_{l=1}^{2h-1}(J-l) [lambda^{2h}] B^J H(J,P;lambda),   [lambda^{2n}]H = N_n(J,P)/prod_{l=1}^{2n-1}(J-l),
# with N_n from HFun.pm (built from the slice series), against the engine polynomials engine_h<h>_r2.txt, h = 1..hmax, as
# polynomial identities in (j1, j2, sigma).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly; use SerS; use HFun;
sub R { Math::BigRat->new(shift) }
my $hmax = shift // 5;
my $N = HFun::numerators($FindBin::Bin, $hmax);
my $J = Poly::add(Poly::var('j1'), Poly::var('j2')); my $P = Poly::mul(Poly::var('j1'), Poly::var('j2'));
sub sp { my $a = shift; my $p = {}; for my $i (0..$#$a) { $p = Poly::add($p, Poly::scale(Poly::pow(Poly::var('s'), $i), $a->[$i])) if $a->[$i] != 0; } return $p; }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
my @bp = (SerS::sone()); for my $p (1..$hmax) { push @bp, SerS::smul($bp[-1], SerS::beta()); }
sub BJ { my $k = shift; my $r = {}; for my $p (0..$k) { $r = Poly::add($r, Poly::mul(sp(SerS::pscale(SerS::coeff($bp[$p], $k), R(1)/fact($p))), Poly::pow($J, $p))); } return $r; }
sub JP { my $q = shift; my $r = {}; for my $k (keys %$q) { my $e = Poly::mono_parse($k); my $t = Poly::const($q->{$k}); for my $v (keys %$e) { my $f = $v eq 'J' ? $J : ($v eq 'P' ? $P : Poly::var($v)); $t = Poly::mul($t, Poly::pow($f, $e->{$v})); } $r = Poly::add($r, $t); } return $r; }
for my $h (1..$hmax) { my $X = CPoly::sigma_form(CPoly::load($FindBin::Bin, $h, 2), $h); my $pred = {};
  for my $n (0..$h) { my $L = 2*$n - 1; $L = 0 if $L < 0; my $t = Poly::mul(JP($N->[$n]), BJ($h - $n));
    for my $l ($L+1 .. 2*$h-1) { $t = Poly::mul($t, Poly::sub_($J, Poly::const($l))); } $pred = Poly::add($pred, $t); }
  print "genus $h: H-form reproduces sigma^h Xi_h(j1,j2) identically: ", (Poly::equal($pred, $X) ? "YES" : "NO"), "\n"; }
