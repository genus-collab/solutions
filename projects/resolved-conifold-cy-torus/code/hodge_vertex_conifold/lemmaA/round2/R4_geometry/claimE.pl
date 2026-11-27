# claimE.pl hmax : the identity E(h,j) = 0 (report, Claim E) for all 1 <= h <= hmax and integers 1 <= j <= 2h-1:
#  E(h,j) = (2h+1)(j-2h) Sum_{n<j/2} chi_{h-n} T_n(j) - Sum_{n<j/2} c_{h+1-n} d/dx (x)_(2h+2-2n) |_{x=j-2n} T_n(j)
#           - Sum_{h''=ceil((j+1)/2)}^{h-1} c_{h+1-h''} (j-2h'')_(2h+2-2h'') Sum_{n<j/2} chi_{h''-n} T_n(j)
#           - j Sum_{n<j/2} zeta'_{h-n} T_n(j),
#  T_n(j) = (j-1)_(2n) [t^{2n}] H(t)^j (Route B's Theorem B data; T_0 = 1),  chi_m = -B_{2m}/(2m),  zeta'_m = (2m-1) B_{2m}/(2m),
#  c_k = 2/(2k)!,  H = S exp(-t S'/S), log H = Sum_k (2k-1) B_{2k} t^{2k}/(2k (2k)!).
# Claim E for all (h,j) implies (by the induction of report s.4) the lambda_1-component of the one-leaf Lemma A in ALL genera.
use strict; use warnings; use Math::BigRat; use Math::BigInt;
sub Q { Math::BigRat->new(shift) }
$| = 1; my $hmax = shift // 20;
my @F = (Q(1)); sub fact { my $n = shift; while (@F <= $n) { push @F, $F[-1] * scalar(@F) } return $F[$n] }
my @B = (Q(1)); sub binom { my ($n, $k) = @_; return Q(0) if $k < 0 || $k > $n; return fact($n) / (fact($k) * fact($n - $k)) }
sub bern { my $n = shift; while ($#B < $n) { my $m = scalar @B; my $s = Q(0); for my $k (0..$m-1) { $s += binom($m + 1, $k) * $B[$k] } push @B, -$s / ($m + 1) } return $B[$n] }
sub chi { my $m = shift; return -bern(2*$m) / (2*$m) }
sub zp { my $m = shift; return (2*$m - 1) * bern(2*$m) / (2*$m) }
sub c { my $k = shift; return Q(2) / fact(2*$k) }
sub ff { my ($x, $n) = @_; my $p = Q(1); $p *= ($x - $_) for 0..$n-1; return $p }
sub dff { my ($x, $n) = @_; my $s = Q(0); for my $i (0..$n-1) { my $p = Q(1); for my $l (0..$n-1) { next if $l == $i; $p *= ($x - $l) } $s += $p } return $s }   # d/dx (x)_(n)
# log H coefficients and H^j at integer j via exp
my $K = $hmax;   # need n < j/2 <= h
my @lh = map { Q(0) } 0..$K; for my $k (1..$K) { $lh[$k] = (2*$k - 1) * bern(2*$k) / (2*$k * fact(2*$k)) }   # coefficient of t^{2k}
my %T;
sub Tn { my ($n, $j) = @_; return $T{"$n|$j"} //= do {
    # [t^{2n}] exp(j log H)
    my @e = (Q(1)); for my $m (1..$n) { my $s = Q(0); for my $k (1..$m) { $s += $k * $lh[$k] * $e[$m - $k] } $e[$m] = $s * $j / $m }
    ff(Q($j) - 1, 2*$n) * $e[$n] } }
my ($tot, $bad) = (0, 0);
for my $h (1..$hmax) { my $hb = 0;
  for my $j (1..2*$h - 1) { my $jq = Q($j); my $E = Q(0);
    for (my $n = 0; 2*$n < $j; $n++) { my $t = Tn($n, $j);
      $E += (2*$h + 1) * ($j - 2*$h) * chi($h - $n) * $t;
      $E -= c($h + 1 - $n) * dff($jq - 2*$n, 2*$h + 2 - 2*$n) * $t;
      for (my $hp = int(($j + 2)/2); $hp <= $h - 1; $hp++) { next unless 2*$hp >= $j + 1; $E -= c($h + 1 - $hp) * ff($jq - 2*$hp, 2*$h + 2 - 2*$hp) * chi($hp - $n) * $t; }
      $E -= $jq * zp($h - $n) * $t; }
    $tot++; if ($E != 0) { $bad++; $hb++; print "  E($h,$j) = $E\n" if $hb <= 3; } }
  printf "genus %2d: Claim E at j = 1..%d: %s\n", $h, 2*$h - 1, ($hb ? "$hb FAILURES" : "all zero") if $h <= 12 || $h % 5 == 0 || $hb; }
printf "Claim E: %d (h,j) pairs checked for h <= %d, %d failures\n", $tot, $hmax, $bad;
