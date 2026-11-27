# fitD2.pl h : test the ansatz (D2) for two leaves,
#   sigma^h Xi_h(j1,j2) = Sum_{m>=0} C(j1-1,m) C(j2-1,m) prod_{l=2m+1}^{2h-1} (J-l) g_{h,m}(J),   g_{h,m} polynomial in J,
# which makes Lemma A (r = 2) evident.  Exact least-squares-free test: solve the linear system coefficientwise (each power of sigma)
# and report whether it is consistent.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly;
sub R { Math::BigRat->new(shift) }
my $h = shift; my $extra = shift // 0;   # extra: allow g_m to depend on P = j1 j2 up to this degree
my $C = CPoly::load($FindBin::Bin, $h, 2); my $S = CPoly::sigma_form($C, $h);
my $d = 3*$h - 1;
sub binomP { my ($v, $m) = @_; my $p = Poly::const(1); for my $i (0..$m-1) { $p = Poly::mul($p, Poly::sub_(Poly::var($v), Poly::const(1 + $i))); } my $f = R(1); $f *= $_ for 2..$m; return Poly::scale($p, R(1)/$f); }
my $J = Poly::add(Poly::var('j1'), Poly::var('j2')); my $P = Poly::mul(Poly::var('j1'), Poly::var('j2'));
my @basis; my @labels;
for my $m (0..int($d/2)) { my $b = Poly::mul(binomP('j1', $m), binomP('j2', $m)); my $pl = Poly::const(1);
  for my $l (2*$m+1 .. 2*$h-1) { $pl = Poly::mul($pl, Poly::sub_($J, Poly::const($l))); }
  my $base = Poly::mul($b, $pl); my $dbase = Poly::total_degree($base, 'j1', 'j2'); next if $dbase > $d;
  for my $e (0..($d - $dbase)) { for my $q (0..$extra) { next if $e + 2*$q > $d - $dbase;
      push @basis, Poly::mul($base, Poly::mul(Poly::pow($J, $e), Poly::pow($P, $q))); push @labels, "m=$m,J^$e,P^$q"; } } }
# monomial index of symmetric polys: use all monomials j1^a j2^b with a+b <= d
my @mons; for my $a (0..$d) { for my $b (0..$d-$a) { push @mons, "$a,$b"; } }
sub coeffvec { my $p = shift; my %c; for my $k (keys %$p) { my $e = Poly::mono_parse($k); my $a = $e->{j1} // 0; my $b = $e->{j2} // 0; $c{"$a,$b"} = $p->{$k}; } return [ map { $c{$_} // R(0) } @mons ]; }
my @A = map { coeffvec($_) } @basis;   # columns
print "genus $h: ansatz with ", scalar(@basis), " basis functions (extra P-degree $extra)\n";
for my $sp (0..$h) { my $target = coeffvec(Poly::coeff_of($S, 's', $sp));
  # solve sum_c x_c A_c = target  (Gaussian elimination on the augmented matrix rows = monomials)
  my $nr = scalar @mons; my $nc = scalar @A; my @M = map { my $i = $_; [ (map { $A[$_][$i] } 0..$nc-1), $target->[$i] ] } 0..$nr-1;
  my $row = 0; my @piv;
  for my $c (0..$nc-1) { my $p; for my $i ($row..$nr-1) { if ($M[$i][$c] != 0) { $p = $i; last } } next unless defined $p;
    @M[$row, $p] = @M[$p, $row]; my $pv = $M[$row][$c]; $_ = $_ / $pv for @{$M[$row]};
    for my $i (0..$nr-1) { next if $i == $row || $M[$i][$c] == 0; my $f = $M[$i][$c]; $M[$i][$_] -= $f * $M[$row][$_] for 0..$nc; }
    push @piv, $c; $row++; }
  my $consistent = 1; for my $i ($row..$nr-1) { $consistent = 0 if $M[$i][$nc] != 0; }
  print "  sigma^$sp: rank $row, ", ($consistent ? "CONSISTENT" : "inconsistent"), "\n";
  if ($consistent) { for my $k (0..$#piv) { my $x = $M[$k][$nc]; print "      $labels[$piv[$k]] : $x\n" if $x != 0; } }
}
