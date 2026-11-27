# Psi.pl nmax : Psi = H/(S1 S2) for two leaves (numerators M_n(J,P) over prod_{l=1}^{2n-1}(J-l)); boundary checks
#   Psi(J, P=0) = 1 (dilaton),  Psi(J, P=J-1) = 1/B (string);  prints M_n by powers of P.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use SerS; use HFun;
sub R { Math::BigRat->new(shift) }
my $nmax = shift // 3;
my $N = HFun::numerators($FindBin::Bin, $nmax);
my $inv = SerS::sexp(SerS::sscale(SerS::logS1S2(), R(-1)));   # 1/(S1S2)
sub sp { my $a = shift; my $p = {}; for my $i (0..$#$a) { $p = Poly::add($p, Poly::scale(Poly::pow(Poly::var('s'), $i), $a->[$i])) if $a->[$i] != 0; } return $p; }
my $J = Poly::var('J');
for my $n (0..$nmax) { my $L = 2*$n - 1; $L = 0 if $L < 0; my $M = {};
  for my $a (0..$n) { my $t = Poly::mul(sp(SerS::coeff($inv, $n - $a)), $N->[$a]); my $La = 2*$a - 1; $La = 0 if $La < 0;
    for my $l ($La+1 .. $L) { $t = Poly::mul($t, Poly::sub_($J, Poly::const($l))); } $M = Poly::add($M, $t); }
  my $M0 = Poly::subst($M, { P => R(0) });
  # string check: P = J - 1
  my $Ms = {}; for my $k (keys %$M) { my $e = Poly::mono_parse($k); my $t = Poly::const($M->{$k}); for my $v (keys %$e) { my $f = $v eq 'P' ? Poly::sub_($J, Poly::const(1)) : Poly::var($v); $t = Poly::mul($t, Poly::pow($f, $e->{$v})); } $Ms = Poly::add($Ms, $t); }
  my $Binv = sp(SerS::coeff(SerS::sexp(SerS::sscale(SerS::beta(), R(-1))), $n)); my $pr = Poly::const(1); for my $l (1..$L) { $pr = Poly::mul($pr, Poly::sub_($J, Poly::const($l))); }
  my $okS = Poly::equal($Ms, $n == 0 ? Poly::const(1) : Poly::mul($Binv, $pr)); my $okD = Poly::equal($M0, $n == 0 ? Poly::const(1) : {});
  print "[lambda^", 2*$n, "] Psi: numerator over prod_{l=1}^{$L}(J-l); dilaton check ", ($okD ? "ok" : "FAIL"), ", string check ", ($okS ? "ok" : "FAIL"), "\n";
  for my $k (0..$n) { my $c = Poly::coeff_of($M, 'P', $k); next if Poly::is_zero($c); print "    P^$k: ", Poly::to_string($c), "\n"; } }
