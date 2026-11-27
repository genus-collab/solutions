# check_vertex_orders.pl -- nontrivial test of the top-relation structure (report.md s.5) against the project's DIRECT genus-2 computation
# (leaf_structure.md, Step 2: bare genus-2 vertex graph, coefficients of t^0 and t^1, exact in x):
#   t^0: -(7 + 40x + 190x^2 - 280x^3 - 184x^4 - 280x^5 + 190x^6 + 40x^7 + 7x^8)/(69120 x c^4)
#   t^1: (1-x)^2 x (5 - 61x + 5x^2)/(3456 c^4)
# Prediction (sub-top cancellation by the string flow; Xi_2(2) = 0; Xi(1..) = 0):
#   [t^{2-2h+k}] = (-x c^3)^{h-1} sum_{D(R) = k} prod_i (-1)^{j_i} C_{j_i} / |Aut R| * Xi^c_h(R),   k = 2, 3 (h = 2),
#   C_j = (c^2/j) eps_j(sigma), Xi^c_h(R) = c^-(3h-3+2r) Xi_h(R)(tau = -x/c).
use strict; use warnings; use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }
for my $x (Q(2), Q(1,3), Q(-3,5)) {
  my $c = 1 + $x; my $tau = -$x / $c; my $sig = $tau * (1 + $tau); my $h = 2;
  my $eps = sub { my $j = shift; my $p = Q(1); for my $k (1..$j-1) { $p *= ($j*$j*$sig + $k*($j-$k)); } return $p / fact($j-1)**2; };
  my $C = sub { my $j = shift; return $c**2/$j * $eps->($j); };
  my $Xc = sub { my @R = @_; my $r = scalar @R; return xi($h, $tau, @R) / $c**(3*$h - 3 + 2*$r); };
  my $pre = (-$x * $c**3)**($h-1);
  my $p0 = $pre * ( -$C->(3) * $Xc->(3) + $C->(2)**2/2 * $Xc->(2,2) );
  my $p1 = $pre * ( $C->(4) * $Xc->(4) - $C->(2)*$C->(3) * $Xc->(2,3) + $C->(2)**3/6 * $Xc->(2,2,2) );
  my $d0 = -(7 + 40*$x + 190*$x**2 - 280*$x**3 - 184*$x**4 - 280*$x**5 + 190*$x**6 + 40*$x**7 + 7*$x**8)/(69120*$x*$c**4);
  my $d1 = (1-$x)**2 * $x * (5 - 61*$x + 5*$x**2)/(3456*$c**4);
  printf "x=%s: t^0 predicted %s, direct %s %s | t^1 predicted %s, direct %s %s\n", $x, $p0, $d0, ($p0 == $d0 ? "ok" : "MISMATCH"), $p1, $d1, ($p1 == $d1 ? "ok" : "MISMATCH");
}
