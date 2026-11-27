package r1lib;
# Round 2, route R1: shared exact (Math::BigRat) helpers.
# Conventions (brief s.1, c = 1): framing tau, sigma = tau(1+tau), Pi_j(tau) = prod_{a=1}^{j-1} (j tau + a)/(j-1)!.
# Conifold leaf curve (leaf_structure_proof.md Theorem 11.1, s.1.5): x = conifold weight, c = 1 + x, framing tau = -x/c,
# conifold coordinate T = 1 - e^{-t}.  A_j(t) = (-1)^j C_j t^-1 zeta(t)^j 2F1(jx/c, jx/c; j; T).
use strict; use warnings;
use Math::BigRat;
use Exporter 'import';
our @EXPORT = qw(R fact binom parts key aut Pi_at series_mul series_inv series_exp series_log series_pow series_compose
                 hyp2f1_series T_of_t_series);

sub R { Math::BigRat->new(@_) }
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; return $f; }
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; return fact($n)/(fact($k)*fact($n-$k)); }
# partitions of n into parts <= max (descending arrays)
sub parts { my ($n, $max) = @_; return ([]) if $n == 0; my @o; for my $k (reverse 1..($max < $n ? $max : $n)) { push @o, [$k, @$_] for parts($n - $k, $k); } return @o; }
sub key { return join(',', sort { $b <=> $a } @_); }
sub aut { my %m; $m{$_}++ for @_; my $a = R(1); $a *= fact($_) for values %m; return $a; }
# Pi_j at a given framing tau (BigRat), returns a new object
sub Pi_at { my ($j, $tau) = @_; my $v = R(1); return $v if $j <= 1; for my $b (1..$j-1) { $v = $v * ($tau*$j + $b); } for my $b (1..$j-1) { $v = $v / $b; } return $v; }

# truncated power series: array refs of BigRat, index = power, length N+1
sub series_mul { my ($a, $b, $N) = @_; my @c = map { R(0) } 0..$N;
  for my $i (0..$N) { next if $i > $#$a || $a->[$i]->is_zero; for my $j (0..$N-$i) { next if $j > $#$b; $c[$i+$j] = $c[$i+$j] + $a->[$i]*$b->[$j]; } } return \@c; }
sub series_inv { my ($a, $N) = @_; die "series_inv: zero constant" if $a->[0]->is_zero; my @b = (R(1)/$a->[0]);
  for my $n (1..$N) { my $s = R(0); for my $k (1..$n) { next if $k > $#$a; $s = $s + $a->[$k]*$b[$n-$k]; } $b[$n] = -$s/$a->[0]; } return \@b; }
sub series_exp { my ($a, $N) = @_; die "series_exp: nonzero constant" unless $a->[0]->is_zero; my @b = (R(1));
  for my $n (1..$N) { my $s = R(0); for my $k (1..$n) { next if $k > $#$a; $s = $s + $k*$a->[$k]*$b[$n-$k]; } $b[$n] = $s/$n; } return \@b; }
sub series_log { my ($a, $N) = @_; die "series_log: constant != 1" unless $a->[0] == 1; my @b = (R(0));
  # b' = a'/a
  my $ia = series_inv($a, $N); my @da = map { ($_+1 <= $#$a) ? ($_+1)*$a->[$_+1] : R(0) } 0..$N; my $q = series_mul(\@da, $ia, $N);
  for my $n (1..$N) { $b[$n] = $q->[$n-1]/$n; } return \@b; }
sub series_pow { my ($a, $e, $N) = @_;   # a[0] = 1, any rational exponent
  my $l = series_log($a, $N); my @m = map { $_ * $e } @$l; return series_exp(\@m, $N); }
# compose f(g(t)) with g(0) = 0
sub series_compose { my ($f, $g, $N) = @_; die "compose: g(0) != 0" unless $g->[0]->is_zero; my @res = map { R(0) } 0..$N; my $pw = [ R(1), map { R(0) } 1..$N ];
  for my $k (0..$N) { last if $k > $#$f; for my $i (0..$N) { $res[$i] = $res[$i] + $f->[$k]*$pw->[$i]; } $pw = series_mul($pw, $g, $N); } return \@res; }
# 2F1(a, b; c; T) as a series in T to order N (a, b, c BigRat; c not a non-positive integer)
sub hyp2f1_series { my ($a, $b, $c, $N) = @_; my @s = (R(1)); my $t = R(1);
  for my $n (1..$N) { $t = $t * ($a + $n - 1) * ($b + $n - 1) / (($c + $n - 1) * $n); push @s, $t; } return \@s; }
# T = 1 - e^{-t} as a series in t
sub T_of_t_series { my $N = shift; my @s = (R(0)); for my $n (1..$N) { push @s, R((-1)**($n+1)) / fact($n); } return \@s; }
1;
