# psi_poly.pl -- the polynomial Psi_R(s) = sum_m 1/m! < Omega prod lam_{j_i} prod_{k<=m} s T(psi_k) >_h  (= Q(-s) of the brief),
# degree <= N0 = 3h-3+r; Xi_h(R) = Psi_R(1). Interpolates Psi_R at s = 0..N0 and reports the multiplicity of the root s = 1.
use strict; use warnings; use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }
sub interp { my ($xs, $ys) = @_; my $n = scalar @$xs; my @c = (Q(0)) x $n;   # Newton -> monomial coefficients
  my @dd = @$ys; my @coef; for my $j (0..$n-1) { push @coef, $dd[$j]; for my $i (reverse $j+1..$n-1) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($xs->[$i] - $xs->[$i-$j-1]); } }
  # rebuild: p(s) = sum coef_j prod_{i<j} (s - x_i)
  my @p = (Q(0)) x $n; my @basis = (Q(1)); for my $j (0..$n-1) { for my $k (0..$#basis) { $p[$k] += $coef[$j] * $basis[$k]; }
    my @nb = (Q(0)) x (scalar(@basis) + 1); for my $k (0..$#basis) { $nb[$k+1] += $basis[$k]; $nb[$k] -= $xs->[$j] * $basis[$k]; } @basis = @nb; }
  return @p; }
sub evalp { my ($p, $x) = @_; my $v = Q(0); for my $k (reverse 0..$#$p) { $v = $v * $x + $p->[$k]; } return $v; }
sub deriv { my @p = @_; my @d; for my $k (1..$#p) { push @d, $p[$k] * $k; } return @d; }
my $tau = Q(2,7);
for my $c ([2,[2]],[2,[3]],[2,[2,2]],[2,[1]],[2,[]],[3,[2]],[3,[3]],[3,[4]],[3,[5]],[3,[2,2]],[3,[2,3]]) {
  my ($h, $R) = @$c; my $r = scalar @$R; my $N0 = 3*$h - 3 + $r;
  my $cls = omega_class($h, Q(1), $tau, -1 - $tau);
  my @xs = map { Q($_) } 0..$N0; my @ys = map { contract($cls, q_comp($h, -$_, @$R)) } @xs;
  my @p = interp(\@xs, \@ys);
  my $mult = 0; my @q = @p; while (@q && evalp(\@q, Q(1)) == 0) { $mult++; @q = deriv(@q); last if $mult > 20; }
  my $D = 0; $D += $_ - 1 for @$R;
  printf "h=%d R=(%s) D=%d N0=%d: deg Psi=%d, mult of root s=1: %d, 2h-2-D=%d, Psi(1)=%s\n", $h, join(",",@$R), $D, $N0, $#p, $mult, 2*$h-2-$D, evalp(\@p, Q(1));
}
