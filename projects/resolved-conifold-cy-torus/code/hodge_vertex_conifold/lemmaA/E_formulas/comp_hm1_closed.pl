# comp_hm1_closed.pl hmax : conjectural closed formula for the component next to lambda_g, two leaves:
#   (-1)^{h-1} Xi^(h-1)_h(j1,j2) = [sigma^1](sigma^h Xi_h)
#     = prod_{l=1}^{2h-1}(J-l) [lambda^{2h}] ( lambda S'(lambda) - J S(lambda) log S(lambda) )
#     + Sum_{k=1}^{h} (j1)_k (j2)_k prod_{l=2k}^{2h-1}(J-l) [lambda^{2h}] S(lambda) (-1)^{k+1} (1/k) Sum_n C(2n-k, k-1) c_n lambda^{2n},
#   S = (l/2)/sinh(l/2), log S = Sum_n c_n l^{2n}, c_n = -B_{2n}/(2n (2n)!).
# (Read off from the sigma-linear parts of the Newton-regular coefficients A_{k,2k-1}, n <= 5; tested here against the engine
#  polynomials engine_h<h>_r2.txt, including genus 6 which was not used to find it.)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly;
sub R { Math::BigRat->new(shift) }
my $hmax = shift // 6;
my @BN = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"));
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $x = R(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } return $x; }
my $K = $hmax;
my @c = (R(0)); for my $n (1..$K) { push @c, -$BN[$n] / (2*$n * fact(2*$n)); }
# S = exp(log S): series coefficients (index n = power of l^2)
my @S = (R(1)); for my $n (1..$K) { my $s = R(0); for my $m (1..$n) { $s += $m * $c[$m] * $S[$n-$m]; } $S[$n] = $s / $n; }
sub smul { my ($a, $b) = @_; my @r = (R(0)) x ($K+1); for my $i (0..$K) { for my $j (0..$K-$i) { $r[$i+$j] += ($a->[$i] // 0) * ($b->[$j] // 0); } } return \@r; }
my @lSp = map { 2*$_ * $S[$_] } 0..$K;                  # l S'(l): coefficient of l^{2n} is 2n S_n
my $SlogS = smul(\@S, \@c);
my $J = Poly::add(Poly::var('j1'), Poly::var('j2'));
sub ff { my ($v, $k) = @_; my $p = Poly::const(1); for my $i (0..$k-1) { $p = Poly::mul($p, Poly::sub_(Poly::var($v), Poly::const($i))); } return $p; }
for my $h (2..$hmax) { my $C = CPoly::load($FindBin::Bin, $h, 2);
  my $data = Poly::scale($C->[$h-1], (($h - 1) % 2) ? R(-1) : R(1));      # (-1)^{h-1} Xi^(h-1)
  my $pred = {}; my $pr = Poly::const(1); $pr = Poly::mul($pr, Poly::sub_($J, Poly::const($_))) for 1..(2*$h-1);
  $pred = Poly::add($pred, Poly::mul($pr, Poly::sub_(Poly::const($lSp[$h]), Poly::scale($J, $SlogS->[$h]))));
  for my $k (1..$h) { my @g = (R(0)) x ($K+1); for my $n ($k..$K) { $g[$n] = ((($k + 1) % 2) ? R(-1) : R(1)) / $k * binom(2*$n - $k, $k - 1) * $c[$n]; }
    my $co = smul(\@S, \@g)->[$h]; next if $co == 0;
    my $t = Poly::mul(ff('j1', $k), ff('j2', $k)); for my $l (2*$k .. 2*$h-1) { $t = Poly::mul($t, Poly::sub_($J, Poly::const($l))); }
    $pred = Poly::add($pred, Poly::scale($t, $co)); }
  print "genus $h: closed formula for Xi^(h-1)_h(j1,j2) ", (Poly::equal($pred, $data) ? "AGREES" : "DISAGREES"), " with the engine polynomial\n"; }
