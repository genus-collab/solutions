# check_lagrange.pl -- exact check of the tree-function (Lagrange) form of Proposition D2.3:
#   Xi_h(R) = (-1)^{N0} [y^{N0}] { (1 - T(y))^{a+N0-1} e^{-(J+N0) T(y)} Phi_R(y) },  T(y) = sum n^{n-1} y^n / n!  (Cayley),
#   Phi_R(y) = sum_p y^p MV_h(R u 1^p)/p!,  a = 2h-2+r, N0 = 3h-3+r.
use strict; use warnings; use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }
sub smul { my ($A, $B, $N) = @_; my @C = (Q(0)) x ($N+1); for my $i (0..$N) { next if $A->[$i] == 0; for my $j (0..$N-$i) { $C[$i+$j] += $A->[$i] * $B->[$j]; } } return \@C; }
sub sexp { my ($A, $N) = @_; # exp of series with A[0] = 0
  my @E = (Q(1), (Q(0)) x $N); my $term = [Q(1), (Q(0)) x $N]; for my $k (1..$N) { $term = smul($term, $A, $N); $term = [ map { $_ / $k } @$term ]; $E[$_] += $term->[$_] for 0..$N; } return \@E; }
sub spow1m { my ($A, $e, $N) = @_; # (1 - A)^e, A[0] = 0, e integer >= 0
  my $R = [Q(1), (Q(0)) x $N]; my $B = [ map { -$_ } @$A ]; $B->[0] = Q(1); for (1..$e) { $R = smul($R, $B, $N); } return $R; }
my $tau = Q(2,7);
for my $c ([2,[2]],[2,[3]],[2,[2,2]],[3,[2]],[3,[3]],[3,[5]]) {
  my ($h, $R) = @$c; my $r = scalar @$R; my $J = 0; $J += $_ for @$R; my $a = 2*$h - 2 + $r; my $N0 = 3*$h - 3 + $r;
  my @T = (Q(0)); for my $n (1..$N0) { push @T, Q($n)**($n-1) / fact($n); }
  my $E = sexp([ map { -($J + $N0) * $_ } @T ], $N0);
  my $P = spow1m(\@T, $a + $N0 - 1, $N0);
  my @Phi = map { mv($h, $tau, @$R, (1) x $_) / fact($_) } 0..$N0;
  my $prod = smul(smul($E, $P, $N0), \@Phi, $N0);
  my $val = Q((-1)**$N0) * $prod->[$N0];
  my $x = xi($h, $tau, @$R);
  printf "h=%d R=(%s): Lagrange form = %s, Xi = %s %s\n", $h, join(",", @$R), $val, $x, ($val == $x ? "ok" : "MISMATCH");
}
