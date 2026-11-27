# fitBB.pl [hmax] : test the "binomial x shortened Norlund" ansatz for two leaves,
#   sigma^h Xi_h(j1,j2) = [lambda^{2h}] B^J Sum_{s1,s2 >= 0} d_{s1,s2}(lambda) C(j1-1,s1) C(j2-1,s2) prod_{l=s1+s2+1}^{2h-1} (J - l),
# with universal (h-independent) series d_{s1,s2} = d_{s2,s1}.  At genus h the residual after subtracting the lower-genus terms
# must have degree <= 2h-1 in (j1,j2) (a genuine test); its decomposition in the genus-h basis then fixes [lambda^{2h}] d.
# Lemma A for r = 2 is evident from this form: C(j1-1,s1)C(j2-1,s2) != 0 needs J >= s1+s2+2, where the product vanishes if J <= 2h-1.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly; use SerS;
sub R { Math::BigRat->new(shift) }
$| = 1;
my $hmax = shift // 5;
my $J = Poly::add(Poly::var('j1'), Poly::var('j2'));
sub binomP { my ($v, $m) = @_; my $p = Poly::const(1); for my $i (0..$m-1) { $p = Poly::mul($p, Poly::sub_(Poly::var($v), Poly::const(1 + $i))); } my $f = R(1); $f *= $_ for 2..$m; return Poly::scale($p, R(1)/$f); }
sub sigpoly { my $a = shift; my $p = {}; for my $i (0..$#$a) { next if $a->[$i] == 0; $p = Poly::add($p, Poly::scale(Poly::pow(Poly::var('s'), $i), $a->[$i])); } return $p; }
my @bp = (SerS::sone()); for my $p (1..$hmax) { push @bp, SerS::smul($bp[-1], SerS::beta()); }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
my %BJ; sub BJ { my $k = shift; return $BJ{$k} //= do { my $r = {}; for my $p (0..$k) { $r = Poly::add($r, Poly::mul(sigpoly(SerS::pscale(SerS::coeff($bp[$p], $k), R(1)/fact($p))), Poly::pow($J, $p))); } $r }; }
my %basis;  # basis{h}{"s1,s2"}
sub basis { my ($h, $s1, $s2) = @_; return $basis{$h}{"$s1,$s2"} //= do { my $p = Poly::mul(binomP('j1', $s1), binomP('j2', $s2)); for my $l ($s1+$s2+1 .. 2*$h-1) { $p = Poly::mul($p, Poly::sub_($J, Poly::const($l))); } $p }; }
my %d;      # d{"s1,s2"}{n} = Poly in s
$d{"0,0"}{0} = Poly::const(1);
for my $h (1..$hmax) {
  my $X = CPoly::sigma_form(CPoly::load($FindBin::Bin, $h, 2), $h);
  for my $key (keys %d) { my ($s1, $s2) = split /,/, $key; for my $n (keys %{$d{$key}}) { next if $n >= $h;
      $X = Poly::sub_($X, Poly::mul(Poly::mul($d{$key}{$n}, basis($h, $s1, $s2)), BJ($h - $n))); } }
  my $deg = Poly::total_degree($X, 'j1', 'j2');
  print "genus $h: residual degree in (j1,j2) = $deg (ansatz needs <= ", 2*$h-1, "): ", ($deg <= 2*$h - 1 ? "OK" : "FAILS"), "\n";
  last if $deg > 2*$h - 1;
  # decompose residual in the basis {C(j1-1,s1)C(j2-1,s2) prod_{l=s1+s2+1}^{2h-1}(J-l) : s1+s2 <= 2h-1}: triangular by decreasing s1+s2 (top binomial part)
  my @keys; for my $t (0..2*$h-1) { for my $s1 (0..$t) { push @keys, [$s1, $t - $s1]; } }
  # linear algebra over the coefficients (monomials in j1,j2 and s)
  my @mons; my %mi; my @cols;
  for my $k (@keys) { my $b = basis($h, @$k); push @cols, $b; for my $m (keys %$b) { $mi{$m} //= scalar(@mons); push @mons, $m unless defined $mons[$mi{$m}]; } }
  for my $m (keys %$X) { my $e = Poly::mono_parse($m); my $sp = $e->{s} // 0; delete $e->{s}; my $mk = Poly::mono_key($e); unless (exists $mi{$mk}) { $mi{$mk} = scalar(@mons); push @mons, $mk; } }
  my $smax = $h + 1; my $ok = 1;
  for my $sp (0..$smax) { my $tgt = Poly::coeff_of($X, 's', $sp);
    my @M = map { my $m = $_; [ (map { $_->{$m} // R(0) } @cols), ($tgt->{$m} // R(0)) ] } @mons;
    my $nc = scalar @cols; my $row = 0; my @piv;
    for my $c (0..$nc-1) { my $p; for my $i ($row..$#M) { if ($M[$i][$c] != 0) { $p = $i; last } } next unless defined $p; @M[$row, $p] = @M[$p, $row];
      my $pv = $M[$row][$c]; $_ = $_ / $pv for @{$M[$row]}; for my $i (0..$#M) { next if $i == $row || $M[$i][$c] == 0; my $f = $M[$i][$c]; $M[$i][$_] -= $f * $M[$row][$_] for 0..$nc; } push @piv, $c; $row++; }
    for my $i ($row..$#M) { $ok = 0 if $M[$i][$nc] != 0; }
    for my $k (0..$#piv) { my ($s1, $s2) = @{$keys[$piv[$k]]}; my $x = $M[$k][$nc]; next if $x == 0;
      $d{"$s1,$s2"}{$h} = Poly::add($d{"$s1,$s2"}{$h} // {}, Poly::scale(Poly::pow(Poly::var('s'), $sp), $x)); } }
  print "   decomposition ", ($ok ? "exact" : "INCONSISTENT"), "\n";
}
for my $key (sort keys %d) { for my $n (sort { $a <=> $b } keys %{$d{$key}}) { my $p = $d{$key}{$n}; next if Poly::is_zero($p); print "  d_{$key} [lambda^", 2*$n, "] = ", join(' + ', map { "(" . Poly::const_val(Poly::coeff_of($p, 's', $_)) . ")s^$_" } grep { Poly::const_val(Poly::coeff_of($p, 's', $_)) != 0 } 0..8), "\n"; } }
