# slice_fit.pl m M : for the slice R = (m, j) of the two-leaf function, test the ansatz
#   sigma^h Xi_h(m, j) = prod_{l=0}^{2h-1-m} (j - l) * [lambda^{2h}] B^j * Sum_{e=0}^{M} j^e C_e(lambda)
# with universal series C_e (independent of h), determined genus by genus (h = 1..6) from the data: at genus h the residual
# after subtracting the lower-genus terms must be a polynomial in j of degree <= M.  (r = 1 analogue: (6.1) with M = 0, C_0 = S1S2.)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use SerS;
sub R { Math::BigRat->new(shift) }
my ($m, $M, $hmax) = @ARGV; $m //= 2; $M //= 1; $hmax //= 5;
# sigma-poly helpers from SerS
my $B = SerS::sexp(SerS::beta());
# [lambda^{2k}] B^j as polynomial in j with sigma-poly coefficients: B^j = exp(j beta): sum_p j^p beta^p/p!
my @bp = (SerS::sone()); for my $p (1..$hmax) { push @bp, SerS::smul($bp[-1], SerS::beta()); }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
sub Bj_coeff { my ($k) = @_; # returns arrayref over p of sigma-polys: [lambda^{2k}] beta^p / p!
  return [ map { SerS::pscale(SerS::coeff($bp[$_], $k), R(1)/fact($_)) } 0..$k ]; }
# data: Q_h(j) evaluated at j = jlist via point evaluation, divided by prod (j-l); fit as polynomial in j
my %C;   # C{e}{n} = sigma-poly
# lambda^0 coefficients from the lambda_g limit: prod_{i=1}^{m-1} (j + i)
{ my @p = (R(1)); for my $i (1..$m-1) { my @q = (R(0)) x (scalar(@p) + 1); for my $k (0..$#p) { $q[$k+1] += $p[$k]; $q[$k] += $p[$k] * $i; } @p = @q; }
  for my $e (0..$M) { $C{$e}{0} = [ $p[$e] // R(0) ]; } }
for my $h (1..$hmax) { Hodge::setup($h); %HodgePt::memo = ();
  my $deg = 3*$h - 3 + 2;  my @js = (2*$h .. 2*$h + $deg + 2);   # points beyond the zeros
  my @vals; for my $j (@js) { my $v = HodgePt::WP([], [$m, $j]); my @s = map { R(0) } 0..$h; for my $k (0..$#$v) { $s[$h - $k] += (($k % 2) ? -1 : 1) * $v->[$k]; }
    my $pr = R(1); $pr *= ($j - $_) for 0..(2*$h - 1 - $m); push @vals, [ map { $_ / $pr } @s ]; }
  # subtract lower-genus contributions: Sum_{n<h} Sum_e C_e,n j^e [lambda^{2(h-n)}] B^j
  my @res = map { [@$_] } @vals;
  for my $n (0..$h-1) { my $bc = Bj_coeff($h - $n);
    for my $e (0..$M) { my $c = $C{$e}{$n} // next;
      for my $ji (0..$#js) { my $j = R($js[$ji]); my $bj = [R(0)]; for my $p (0..$#$bc) { $bj = SerS::padd($bj, SerS::pscale($bc->[$p], $j ** $p)); }
        my $t = SerS::pmul($c, SerS::pscale($bj, $j ** $e)); for my $sp (0..$#$t) { $res[$ji][$sp] = ($res[$ji][$sp] // R(0)) - $t->[$sp]; } } } }
  # interpolate residual as polynomial in j for each sigma power; check degree <= M
  my $ok = 1; my @coef;
  for my $sp (0..$h + 2) { my @y = map { $_->[$sp] // R(0) } @res;
    # Newton divided differences -> polynomial coefficients
    my @x = map { R($_) } @js; my @dd = @y; my $n = scalar @x;
    for my $k (1..$n-1) { for (my $i = $n-1; $i >= $k; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($x[$i] - $x[$i-$k]); } }
    my @p = (R(0)) x $n; my @basis = (R(1)); # build polynomial from Newton form
    for my $k (0..$n-1) { for my $i (0..$#basis) { $p[$i] += $dd[$k] * $basis[$i]; } my @nb = (R(0)) x (scalar(@basis) + 1); for my $i (0..$#basis) { $nb[$i+1] += $basis[$i]; $nb[$i] -= $basis[$i] * $x[$k]; } @basis = @nb; }
    for my $i ($M+1 .. $#p) { $ok = 0 if $p[$i] != 0; }
    for my $e (0..$M) { $coef[$e][$sp] = $p[$e] // R(0); } }
  for my $e (0..$M) { $C{$e}{$h} = $coef[$e]; }
  print "genus $h (slice m = $m, M = $M): residual ", ($ok ? "has degree <= $M: consistent" : "has degree > $M: ANSATZ FAILS"), "\n";
  last unless $ok;
}
for my $e (0..$M) { for my $n (sort { $a <=> $b } keys %{$C{$e}}) { my $p = $C{$e}{$n}; print "  C_$e [lambda^", 2*$n, "] = ", join(' ', map { "($p->[$_])s^$_" } grep { $p->[$_] != 0 } 0..$#$p), "\n"; } }
