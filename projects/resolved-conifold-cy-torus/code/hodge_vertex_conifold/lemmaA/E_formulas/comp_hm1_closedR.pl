# comp_hm1_closedR.pl r hmax : test the pairwise extension of the lambda_{h-1}-component formula to r leaves:
#   [sigma^1](sigma^h Xi_h(R)) = prod_{l=1}^{N}(J-l)[l^{2h}](l S' - J S log S)
#        + Sum_{i<i'} Sum_{k>=1} (j_i)_k (j_i')_k prod_{l=2k}^{N}(J-l) [l^{2h}] S g_k,   g_k = (-1)^{k+1}(1/k) Sum_n C(2n-k,k-1) c_n l^{2n}.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly;
sub R { Math::BigRat->new(shift) }
my ($r, $hmax) = @ARGV;
my @BN = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"));
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $x = R(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } return $x; }
my $K = $hmax; my @c = (R(0)); for my $n (1..$K) { push @c, -$BN[$n] / (2*$n * fact(2*$n)); }
my @S = (R(1)); for my $n (1..$K) { my $s = R(0); for my $m (1..$n) { $s += $m * $c[$m] * $S[$n-$m]; } $S[$n] = $s / $n; }
sub smul { my ($a, $b) = @_; my @q = (R(0)) x ($K+1); for my $i (0..$K) { for my $j (0..$K-$i) { $q[$i+$j] += ($a->[$i] // 0) * ($b->[$j] // 0); } } return \@q; }
my @lSp = map { 2*$_ * $S[$_] } 0..$K; my $SlogS = smul(\@S, \@c);
my $J = {}; $J = Poly::add($J, Poly::var("j$_")) for 1..$r;
sub ff { my ($v, $k) = @_; my $p = Poly::const(1); for my $i (0..$k-1) { $p = Poly::mul($p, Poly::sub_(Poly::var($v), Poly::const($i))); } return $p; }
for my $h (2..$hmax) { my $C = CPoly::load($FindBin::Bin, $h, $r); my $N = 2*$h - 3 + $r;
  my $data = Poly::scale($C->[$h-1], (($h - 1) % 2) ? R(-1) : R(1));
  my $pr = Poly::const(1); $pr = Poly::mul($pr, Poly::sub_($J, Poly::const($_))) for 1..$N;
  my $pred = Poly::mul($pr, Poly::sub_(Poly::const($lSp[$h]), Poly::scale($J, $SlogS->[$h])));
  for my $k (1..$h) { my @g = (R(0)) x ($K+1); for my $n ($k..$K) { $g[$n] = ((($k + 1) % 2) ? R(-1) : R(1)) / $k * binom(2*$n - $k, $k - 1) * $c[$n]; }
    my $co = smul(\@S, \@g)->[$h]; next if $co == 0;
    my $pp = Poly::const(1); $pp = Poly::mul($pp, Poly::sub_($J, Poly::const($_))) for (2*$k)..$N;
    for my $i (1..$r) { for my $i2 ($i+1..$r) { $pred = Poly::add($pred, Poly::scale(Poly::mul(Poly::mul(ff("j$i", $k), ff("j$i2", $k)), $pp), $co)); } } }
  my $diff = Poly::sub_($data, $pred);
  print "genus $h, r = $r: pairwise formula ", (Poly::is_zero($diff) ? "AGREES" : "DISAGREES (residual total degree " . Poly::total_degree($diff) . ")"), "\n"; }
