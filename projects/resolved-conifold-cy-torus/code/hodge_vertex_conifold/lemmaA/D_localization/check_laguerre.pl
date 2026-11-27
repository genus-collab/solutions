# check_laguerre.pl -- exact check of the dilaton-leaf / degree-1-edge identity (report.md, Proposition D2.3):
#   Xi_h(R) = sum_{p=0}^{N0} (-1)^p L^{(a+p)}_{N0-p}(-(J+p)) MV_h(R u 1^p) / p!,   a = 2h-2+r, N0 = 3h-3+r,
#   MV_h(mu) = Int_{Mbar_{h,l}} Omega prod lam_{mu_i} (no dilaton leaves), L^{(alpha)}_n the Laguerre polynomial.
use strict; use warnings; use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }
sub laguerre { my ($n, $al, $x) = @_; my $t = Q(0); for my $k (0..$n) { $t += binom($n + $al, $n - $k) * (-$x)**$k / fact($k); } return $t; }
my @cases = ([2, [2]], [2, [3]], [2, [2,2]], [2, [1]], [3, [2]], [3, [3]], [3, [4]], [3, [2,2]], [3, [5]]);
for my $c (@cases) { my ($h, $R) = @$c; my $r = scalar @$R; my $J = 0; $J += $_ for @$R;
  my $a = 2*$h - 2 + $r; my $N0 = 3*$h - 3 + $r;
  for my $tau (Q(2,7), Q(-3,5)) {
    my $sum = Q(0);
    for my $p (0..$N0) { my $phi = mv($h, $tau, @$R, (1) x $p); $sum += Q((-1)**$p) * laguerre($N0 - $p, $a + $p, Q(-($J + $p))) * $phi / fact($p); }
    my $x = xi($h, $tau, @$R);
    printf "h=%d R=(%s) tau=%s : Laguerre sum = %s , Xi = %s  %s\n", $h, join(",", @$R), $tau, $sum, $x, ($sum == $x ? "ok" : "MISMATCH");
  } }
