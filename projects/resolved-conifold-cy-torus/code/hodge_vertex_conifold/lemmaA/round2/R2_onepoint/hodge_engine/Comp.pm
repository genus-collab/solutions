package Comp;
# Conversion of Hodge-side data (c-normalisation, symbol x) to the c = 1 normalisation and to the tau-independent
# components:  Xi^{c=1}(tau) = c^{3h-3+2r} Xi^{c}(x),  c = 1+x, tau = -x/(1+x), sigma = tau(1+tau) = -x/c^2, u = 1/sigma,
#              Xi^{c=1} = Sum_{k=0}^{h} (-u)^k Xi^{(k)}   (proof_status_mixed_scale.md s.1.1(e)).
use strict; use warnings;
use Math::BigRat; use Poly; use MParse;
sub R { Math::BigRat->new(shift) }
# $expr: Mathematica string in x and other symbols; returns arrayref of component polys [Xi^(0), ..., Xi^(h)]
# and a consistency flag (checked at extra points).
sub components {
  my ($expr, $h, $r, %opt) = @_;
  my @xs = map { R($_) } ("1/2", "3", "5", "2/7", "-3", "-1/5", "11", "4", "1/9", "13/4");
  my $np = $h + 1 + ($opt{extra} // 2);
  my (@us, @vals);
  for my $i (0..$np-1) { my $x = $xs[$i]; my $c = 1 + $x;
    my $p = MParse::parse($expr, { x => $x });
    $p = Poly::scale($p, $c ** (3*$h - 3 + 2*$r));
    push @us, -($c*$c)/$x; push @vals, $p; }
  # interpolate in v = -u : Xi = Sum_k v^k Xi^(k), using the first h+1 points
  my @vs = map { -$_ } @us;
  my @comp = interp([@vs[0..$h]], [@vals[0..$h]], $h);
  my $ok = 1;
  for my $i ($h+1 .. $np-1) { my $pred = {}; for my $k (0..$h) { $pred = Poly::add($pred, Poly::scale($comp[$k], $vs[$i] ** $k)); }
    $ok = 0 unless Poly::equal($pred, $vals[$i]); }
  return (\@comp, $ok);
}
# Newton/Lagrange: given points v_i and polynomial values P_i, find coefficients C_k (k = 0..d) with Sum_k v_i^k C_k = P_i
sub interp { my ($v, $P, $d) = @_; my $n = $d + 1;
  # solve Vandermonde by Lagrange basis expansion
  my @C = map { {} } 0..$d;
  for my $i (0..$d) { # basis poly L_i(v) = prod_{m != i} (v - v_m)/(v_i - v_m), as coefficient list
    my @L = (R(1)); my $den = R(1);
    for my $m (0..$d) { next if $m == $i; my @N = map { R(0) } 0..$#L+1;
      for my $k (0..$#L) { $N[$k+1] += $L[$k]; $N[$k] -= $L[$k] * $v->[$m]; } @L = @N; $den *= ($v->[$i] - $v->[$m]); }
    for my $k (0..$d) { $C[$k] = Poly::add($C[$k], Poly::scale($P->[$i], $L[$k] / $den)) if $L[$k] != 0; } }
  return @C; }
1;
