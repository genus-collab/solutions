# comp_hm1_r3fit.pl hmax : sigma-linear part of the three-leaf function minus the J-term and the pairwise terms;
# fit the residual in the basis m_k(j) prod_{l=|k|}^{N}(J-l), k = (k1>=k2>=k3>=1) (triple contractions), report coefficients.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly;
sub R { Math::BigRat->new(shift) }
my $hmax = shift // 5; my $r = 3;
my @BN = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"));
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $x = R(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } return $x; }
my $K = $hmax; my @c = (R(0)); for my $n (1..$K) { push @c, -$BN[$n] / (2*$n * fact(2*$n)); }
my @S = (R(1)); for my $n (1..$K) { my $s = R(0); for my $m (1..$n) { $s += $m * $c[$m] * $S[$n-$m]; } $S[$n] = $s / $n; }
sub smul { my ($a, $b) = @_; my @q = (R(0)) x ($K+1); for my $i (0..$K) { for my $j (0..$K-$i) { $q[$i+$j] += ($a->[$i] // 0) * ($b->[$j] // 0); } } return \@q; }
my @lSp = map { 2*$_ * $S[$_] } 0..$K; my $SlogS = smul(\@S, \@c);
my $J = {}; $J = Poly::add($J, Poly::var("j$_")) for 1..$r;
sub ff { my ($v, $k) = @_; my $p = Poly::const(1); for my $i (0..$k-1) { $p = Poly::mul($p, Poly::sub_(Poly::var($v), Poly::const($i))); } return $p; }
sub perms { my @x = @_; return ([]) unless @x; my %seen; my @out; for my $i (0..$#x) { next if $seen{$x[$i]}++; my @rest = @x; splice(@rest, $i, 1); push @out, [$x[$i], @$_] for perms(@rest); } return @out; }
sub mk { my @k = @_; my $p = {}; for my $q (perms(@k)) { my $t = Poly::const(1); for my $i (0..$r-1) { $t = Poly::mul($t, ff("j" . ($i+1), $q->[$i])); } $p = Poly::add($p, $t); } return $p; }
for my $h (2..$hmax) { my $C = CPoly::load($FindBin::Bin, $h, $r); my $N = 2*$h - 3 + $r;
  my $res = Poly::scale($C->[$h-1], (($h - 1) % 2) ? R(-1) : R(1));
  my $pr = Poly::const(1); $pr = Poly::mul($pr, Poly::sub_($J, Poly::const($_))) for 1..$N;
  $res = Poly::sub_($res, Poly::mul($pr, Poly::sub_(Poly::const($lSp[$h]), Poly::scale($J, $SlogS->[$h]))));
  for my $k (1..$h) { my @g = (R(0)) x ($K+1); for my $n ($k..$K) { $g[$n] = ((($k + 1) % 2) ? R(-1) : R(1)) / $k * binom(2*$n - $k, $k - 1) * $c[$n]; }
    my $co = smul(\@S, \@g)->[$h]; next if $co == 0; my $pp = Poly::const(1); $pp = Poly::mul($pp, Poly::sub_($J, Poly::const($_))) for (2*$k)..$N;
    for my $i (1..$r) { for my $i2 ($i+1..$r) { $res = Poly::sub_($res, Poly::scale(Poly::mul(Poly::mul(ff("j$i", $k), ff("j$i2", $k)), $pp), $co)); } } }
  # basis of triple contractions
  my @lab; my @bas;
  for my $a (1..$h+2) { for my $b (1..$a) { for my $cc (1..$b) { my $kk = $a + $b + $cc; next if $kk > $N + 1; my $pp = Poly::const(1); $pp = Poly::mul($pp, Poly::sub_($J, Poly::const($_))) for $kk..$N;
        push @lab, "($a,$b,$cc)"; push @bas, Poly::mul(mk($a, $b, $cc), $pp); } } }
  my %mi; my @mons; for my $b (@bas, $res) { for my $m (keys %$b) { unless (exists $mi{$m}) { $mi{$m} = scalar @mons; push @mons, $m; } } }
  my @rows = map { my $m = $_; [ (map { $_->{$m} // R(0) } @bas), ($res->{$m} // R(0)) ] } @mons;
  my $nc = scalar @bas; my $row = 0; my @piv;
  for my $cidx (0..$nc-1) { my $p; for my $i ($row..$#rows) { if ($rows[$i][$cidx] != 0) { $p = $i; last } } next unless defined $p; @rows[$row, $p] = @rows[$p, $row];
    my $pv = $rows[$row][$cidx]; $_ = $_ / $pv for @{$rows[$row]}; for my $i (0..$#rows) { next if $i == $row || $rows[$i][$cidx] == 0; my $f = $rows[$i][$cidx]; $rows[$i][$_] -= $f * $rows[$row][$_] for 0..$nc; } push @piv, $cidx; $row++; }
  my $ok = 1; for my $i ($row..$#rows) { $ok = 0 if $rows[$i][$nc] != 0; }
  print "genus $h: residual ", (Poly::is_zero($res) ? "zero" : "non-zero"), "; triple-contraction fit ", ($ok ? "CONSISTENT" : "INCONSISTENT"), " (rank $row of $nc)\n";
  if ($ok) { for my $q (0..$#piv) { my $x = $rows[$q][$nc]; print "    $lab[$piv[$q]]: $x\n" if $x != 0; } } }
