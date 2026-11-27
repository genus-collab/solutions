# candK1.pl h : test the kappa-coordinate "shifted lambda_g" candidate
#   sigma^h Int_{Mbar_h} Omega K kappa_{c_1}...kappa_{c_s}  =?  [lambda^{2h}] S1 S2  L0( prod_i (e^{beta d} kappa)_{c_i} ),
#   (e^{beta d} kappa)_c = Sum_{m=0}^{c} beta^m/m! kappa_{c-m},  L0(X) = Int_{Mbar_{0,2h}} K X  (2h points without insertions),
# for s = 1 (this is (6.1)) and s = 2, 3.  Prints the mismatching kappa monomials.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use SerS;
sub R { Math::BigRat->new(shift) }
my $h = $main::h // 2; my $smax = 0;
Hodge::setup($h); my $D = 3*$h - 3;
my $Kx = Hodge::kexp($D + 2*$h);
sub L0 { my @beta = grep { $_ > 0 } @_; my $c0 = R(1); for (@_) { $c0 *= (2*$h - 2) if $_ == 0 } return R(0) if grep { $_ < 0 } @_;
  my $n = 2*$h; my $res = R(0); my $db = 0; $db += $_ for @beta;
  for my $gk (keys %$Kx) { my @gm = $gk eq '' ? () : split /,/, $gk; my $dg = 0; $dg += $_ for @gm; next unless $db + $dg == $n - 3;
    $res += $Kx->{$gk} * Hodge::HK0(0, [@beta, @gm], [(0) x $n]); }
  return $res * $c0; }
my $S12 = SerS::S1S2(); my $bet = SerS::beta();
my @bpow = (SerS::sone()); for my $m (1..$D) { push @bpow, SerS::smul($bpow[-1], $bet); }
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 2..$n; $r }
my @coefM = map { SerS::coeff(SerS::smul($S12, $bpow[$_]), $h) } 0..$D;   # [L^h] S1S2 beta^M (sigma-poly)
sub lhs { my @c = @_; my $v = Hodge::T(@c);   # v-poly, v = -1/sigma ; sigma^h * value = Sum_k (-1)^k sigma^(h-k) T^(k)
  my @p = map { R(0) } 0..$h; for my $k (0..$#$v) { $p[$h - $k] += (($k % 2) ? -1 : 1) * $v->[$k]; } return \@p; }
sub rhs { my @c = @_; my $tot = [R(0)];
  my @ms = ([]); for my $ci (@c) { my @new; for my $m (@ms) { for my $x (0..$ci) { push @new, [@$m, $x]; } } @ms = @new; }
  for my $m (@ms) { my $M = 0; $M += $_ for @$m; my $w = R(1); $w /= fact($_) for @$m;
    my $l0 = L0(map { $c[$_] - $m->[$_] } 0..$#c); next if $l0 == 0;
    $tot = SerS::padd($tot, SerS::pscale($coefM[$M], $w * $l0)); }
  return $tot; }
sub peq { my ($a, $b) = @_; my $n = $#$a > $#$b ? $#$a : $#$b; for (0..$n) { return 0 if ($a->[$_] // 0) != ($b->[$_] // 0) } return 1; }
for my $s (1..$smax) { my ($n, $bad) = (0, 0); my @fails;
  for my $p (Hodge::partitions($D)) { }   # (kappa monomials of all degrees <= D)
  for my $d (0..$D) { for my $p (Hodge::partitions($d)) { next unless scalar(@$p) == $s; $n++;
      my $l = lhs(@$p); my $r = rhs(@$p); unless (peq($l, $r)) { $bad++; push @fails, "{@$p}"; } } }
  # also monomials containing kappa_0 are trivial multiples; skipped
  print "genus $h, s = $s: $n kappa monomials tested, $bad mismatches", ($bad ? " (first: @fails[0..($#fails < 5 ? $#fails : 5)])" : ""), "\n"; }
sub pstr { my $p = shift; my @t; for my $i (0..$#$p) { next if $p->[$i] == 0; push @t, "($p->[$i])*s^$i"; } return @t ? join(' + ', @t) : '0'; }
for my $d (0..$D) { for my $p (Hodge::partitions($d)) { next unless scalar(@$p) == 2;
  my $l = lhs(@$p); my $r = rhs(@$p); my $df = SerS::padd($l, SerS::pscale($r, R(-1)));
  print "{@$p}: lhs = ", pstr($l), "\n       rhs = ", pstr($r), "\n      diff = ", pstr($df), "\n"; } }
