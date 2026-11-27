package Data;
# Reading the exact Hodge-side tables of lemmaA/data and extracting the sigma -> infinity ("Hurwitz", k = 0) component.
# A table value E(x,...) is in the c-normalisation (c = 1+x); for r leaves (probes) the c = 1 value is
#   F = c^(3g-3+2r) E,  and F = sum_{k=0}^{g} t^k F^(k),  t = -1/sigma = c^2/x   (proof_status_mixed_scale.md s.1.1(e)).
# The k = 0 component (single class Lambda^vee(1)) is F at t = 0; we get it by Lagrange interpolation in t from x = 1..g+1
# and check the degree bound with two extra points x = g+2, g+3.
use strict; use warnings; use Math::BigRat; use HB;
sub BR { Math::BigRat->new(shift) }
my $DATADIR = "../data";
sub load { my ($name) = @_; my $recs = HB::read_record("$DATADIR/$name"); my @out;
  for my $r (@$recs) { my $expr = $r->[-1]; my $ast = HB::fold(HB::parse($expr)); push @out, { fields => $r, ast => $ast } }
  return \@out }
sub interp0 { my ($ts, $vs) = @_; my $tot = BR(0);
  for my $i (0..$#$ts) { my $w = $vs->[$i]; for my $k (0..$#$ts) { next if $k == $i; $w *= (BR(0) - $ts->[$k]) / ($ts->[$i] - $ts->[$k]) } $tot += $w }
  return $tot }
sub interp_at { my ($ts, $vs, $t0) = @_; my $tot = BR(0);
  for my $i (0..$#$ts) { my $w = $vs->[$i]; for my $k (0..$#$ts) { next if $k == $i; $w *= ($t0 - $ts->[$k]) / ($ts->[$i] - $ts->[$k]) } $tot += $w }
  return $tot }
# components F^(0..g) of a table value at fixed other variables (env without x); returns arrayref of g+1 BigRats
sub components { my ($ast, $g, $r, $env) = @_; my (@ts, @vs);
  for my $xv (1..$g + 3) { my %e = (%$env, x => BR($xv)); my $E = HB::evaluate($ast, \%e); my $c = BR(1 + $xv);
    push @ts, $c * $c / BR($xv); push @vs, $E * ($c ** (3*$g - 3 + 2*$r)) }
  my @t = @ts[0..$g]; my @v = @vs[0..$g];
  for my $extra ($g + 1, $g + 2) { my $pred = interp_at(\@t, \@v, $ts[$extra]); die "degree check failed in t (not a polynomial of degree <= $g)" unless $pred == $vs[$extra] }
  # coefficients of the polynomial in t: solve via Newton divided differences then expand
  my @coef = poly_from_points(\@t, \@v); return \@coef }
sub poly_from_points { my ($t, $v) = @_; my $n = scalar @$t; my @dd = map { $_->copy } @$v;
  for my $j (1..$n-1) { for (my $i = $n-1; $i >= $j; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($t->[$i] - $t->[$i-$j]) } }
  my @p = (BR(0)) x $n; # Horner expansion of Newton form
  my @acc = (BR(0)) x $n;
  for (my $i = $n-1; $i >= 0; $i--) { # acc = acc*(T - t_i) + dd_i
    my @new = (BR(0)) x $n; for my $k (0..$n-2) { $new[$k+1] += $acc[$k]; } for my $k (0..$n-1) { $new[$k] -= $acc[$k] * $t->[$i] } $new[0] += $dd[$i]; @acc = @new }
  return @acc }
1;
