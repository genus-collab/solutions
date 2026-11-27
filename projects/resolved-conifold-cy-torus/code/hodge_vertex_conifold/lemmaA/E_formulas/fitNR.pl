# fitNR.pl r hmax [diag] : test the Newton-regular (NR) form for r leaves directly on the engine polynomials:
#   sigma^h Xi_h(R) = Sum_{k, n'} [lambda^{2h}] ( S1 S2 B^J A_{k,n'}(lambda) ) m_k(j) prod_{l=n'+1}^{N} (J - l),   N = 2h-3+r,
# m_k = symmetrisation of prod_i (j_i)_{k_i} over the partition k (|k| >= 1), A_{0,0} = 1, and 0 <= n' <= |k| - 1 (this range makes
# every term vanish on the Lemma-A lattice range).  A_{k,n'} universal series (independent of h), determined genus by genus;
# at genus h the residual must lie in the span of the genus-h basis {m_k prod_{l=n'+1}^{N}(J-l)}.
# With diag = 1 (r = 2) only k = (k,k) is allowed (the form found in PsiNewton.pl).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly; use SerS;
sub R { Math::BigRat->new(shift) }
$| = 1;
my ($r, $hmax, $diag) = @ARGV; $diag //= 0;
my @jv = map { "j$_" } 1..$r;
my $J = {}; $J = Poly::add($J, Poly::var($_)) for @jv;
sub sp { my $a = shift; my $p = {}; for my $i (0..$#$a) { $p = Poly::add($p, Poly::scale(Poly::pow(Poly::var('s'), $i), $a->[$i])) if $a->[$i] != 0; } return $p; }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
my @bp = (SerS::sone()); for my $p (1..$hmax) { push @bp, SerS::smul($bp[-1], SerS::beta()); }
my $S12 = SerS::S1S2();
my %SB; sub SB { my $k = shift; return $SB{$k} //= do {   # [lambda^{2k}] S1S2 B^J as Poly in J (j's) and s
    my $r = {}; for my $a (0..$k) { my $c = SerS::coeff($S12, $k - $a); for my $p (0..$a) { my $bc = SerS::pscale(SerS::coeff($bp[$p], $a), R(1)/fact($p));
        $r = Poly::add($r, Poly::mul(sp(SerS::pmul($c, $bc)), Poly::pow($J, $p))); } } $r }; }
sub ffall { my ($v, $k) = @_; my $p = Poly::const(1); for my $i (0..$k-1) { $p = Poly::mul($p, Poly::sub_(Poly::var($v), Poly::const($i))); } return $p; }
sub perms { my @x = @_; return ([]) unless @x; my %seen; my @out; for my $i (0..$#x) { next if $seen{$x[$i]}++; my @rest = @x; splice(@rest, $i, 1); push @out, [$x[$i], @$_] for perms(@rest); } return @out; }
my %mk; sub mk { my @k = @_; my $key = join(',', @k); return $mk{$key} //= do { my $p = {}; for my $q (perms(@k)) { my $t = Poly::const(1); for my $i (0..$r-1) { $t = Poly::mul($t, ffall($jv[$i], $q->[$i])); } $p = Poly::add($p, $t); } $p }; }
sub partitions { my ($n, $parts, $max) = @_; return ([(0) x $parts]) if $n == 0; return () if $parts == 0; my @out; $max //= $n;
  for my $a (reverse 1..($max < $n ? $max : $n)) { for my $q (partitions($n - $a, $parts - 1, $a)) { push @out, [$a, @$q]; } } return @out; }
my %A;   # A{"k|n'"}{n} = Poly in s
$A{"0|0"}{0} = Poly::const(1);
for my $h (1..$hmax) { my $N = 2*$h - 3 + $r;
  my $X = CPoly::sigma_form(CPoly::load($FindBin::Bin, $h, $r), $h);
  for my $key (keys %A) { my ($ks, $np) = split /\|/, $key; my @k = split /,/, $ks; my $b = Poly::mul(mk(@k), do { my $p = Poly::const(1); $p = Poly::mul($p, Poly::sub_($J, Poly::const($_))) for ($np+1)..$N; $p });
    for my $n (keys %{$A{$key}}) { next if $n >= $h; $X = Poly::sub_($X, Poly::mul(Poly::mul($A{$key}{$n}, $b), SB($h - $n))); } }
  # genus-h basis
  my @lab; my @bas;
  for my $K (1 .. $h + $N) { for my $k (partitions($K, $r)) { next if $diag && $r == 2 && $k->[0] != $k->[1];
      for my $np (0 .. $K - 1) { next if $np > $N; my $deg = $K + $N - $np; next if $deg > 3*$h - 3 + $r;
        my $p = Poly::const(1); $p = Poly::mul($p, Poly::sub_($J, Poly::const($_))) for ($np+1)..$N;
        push @lab, join(',', @$k) . "|$np"; push @bas, Poly::mul(mk(@$k), $p); } } }
  # solve per sigma power
  my %mi; my @mons; for my $b (@bas, $X) { for my $m (keys %$b) { my $e = Poly::mono_parse($m); delete $e->{s}; my $kk = Poly::mono_key($e); unless (exists $mi{$kk}) { $mi{$kk} = scalar @mons; push @mons, $kk; } } }
  my $ok = 1; my $rankmax = 0;
  for my $spow (0..$h + 1) { my $tgt = Poly::coeff_of($X, 's', $spow); next if Poly::is_zero($tgt) && $spow > 0;
    my @rows = map { my $m = $_; [ (map { $_->{$m} // R(0) } @bas), ($tgt->{$m} // R(0)) ] } @mons;
    my $nc = scalar @bas; my $row = 0; my @piv;
    for my $c (0..$nc-1) { my $p; for my $i ($row..$#rows) { if ($rows[$i][$c] != 0) { $p = $i; last } } next unless defined $p; @rows[$row, $p] = @rows[$p, $row];
      my $pv = $rows[$row][$c]; $_ = $_ / $pv for @{$rows[$row]}; for my $i (0..$#rows) { next if $i == $row || $rows[$i][$c] == 0; my $f = $rows[$i][$c]; $rows[$i][$_] -= $f * $rows[$row][$_] for 0..$nc; } push @piv, $c; $row++; }
    for my $i ($row..$#rows) { $ok = 0 if $rows[$i][$nc] != 0; }
    $rankmax = $row if $row > $rankmax;
    for my $q (0..$#piv) { my $x = $rows[$q][$nc]; next if $x == 0; $A{$lab[$piv[$q]]}{$h} = Poly::add($A{$lab[$piv[$q]]}{$h} // {}, Poly::scale(Poly::pow(Poly::var('s'), $spow), $x)); } }
  printf "genus %d, r = %d: %d basis functions (rank %d): %s\n", $h, $r, scalar(@bas), $rankmax, ($ok ? "CONSISTENT" : "INCONSISTENT"); last unless $ok; }
