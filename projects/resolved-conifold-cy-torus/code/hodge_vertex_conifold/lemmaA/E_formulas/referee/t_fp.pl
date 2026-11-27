# t_fp.pl gmax : the literature statement cited in report s.4(i), Faber-Pandharipande (Invent. Math. 139 (2000)), Theorem 2:
#   Sum_{g >= 0} t^{2g} Sum_{i=0}^{g} k^i Int_{Mbar_{g,1}} psi^{2g-2+i} lambda_{g-i} = ((t/2)/sin(t/2))^{k+1},
# checked with the referee engine for g <= gmax, k = 0..5 (normalisation check of the citation).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
my $gmax = shift // 5;
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
my $K = $gmax; my @sn = map { ((($_ % 2) ? Q(-1) : Q(1)) / (Q(4) ** $_) / fact(2*$_ + 1)) } 0..$K; my @inv = (Q(1));
for my $n (1..$K) { my $s = Q(0); for my $m (1..$n) { $s -= $sn[$m] * $inv[$n - $m]; } $inv[$n] = $s; }
sub spow { my ($a, $e) = @_; my @r = (Q(1), map { Q(0) } 1..$K); for (1..$e) { my @n = map { Q(0) } 0..$K; for my $i (0..$K) { for my $j (0..$K-$i) { $n[$i+$j] += $r[$i] * $a->[$j]; } } @r = @n; } return \@r; }
my ($n, $bad) = (0, 0);
for my $g (1..$gmax) { my $maxd = $g == 1 ? 1 : 3*$g - 3; my @lam = RefHodge::lambdas($g, $maxd < 1 ? 1 : $maxd);
  my @I; for my $i (0..$g) { my $s = Q(0); for my $key (keys %{$lam[$g - $i]}) { my @mu = split /,/, $key; $s += $lam[$g - $i]{$key} * RefHodge::I($g, \@mu, [2*$g - 2 + $i]); } push @I, $s; }
  for my $k (0..5) { my $lhs = Q(0); $lhs += Q($k) ** $_ * $I[$_] for 0..$g; my $rhs = spow(\@inv, $k + 1)->[$g]; $n++; $bad++ if $lhs != $rhs; } }
printf "Faber-Pandharipande Thm 2 vs referee engine: %d comparisons (g = 1..%d, k = 0..5), %d mismatches\n", $n, $gmax, $bad;
