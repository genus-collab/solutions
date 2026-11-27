# binom_basis.pl h r : binomial-basis coefficients alpha^(k)_a of the components,
#   Xi^(k)_h(j_1..j_r) = Sum_a alpha^(k)_a prod_i C(j_i - 1, a_i),
# computed by Moebius inversion from the values at the lattice points j_i - 1 = b_i >= 0:
#   alpha_a = Sum_{b <= a} (-1)^{|a-b|} prod C(a_i, b_i) Xi(b + 1).
# Lemma A (componentwise)  <=>  alpha^(k)_a = 0 for |a| <= 2h-3.  Prints the table by layers |a| = 2h-2, ..., 3h-3+r.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Comps;
sub R { Math::BigRat->new(shift) }
my ($h, $r) = @ARGV; my $C = Comps::load("$FindBin::Bin/engine_h${h}_r${r}.txt");
my $N = 3*$h - 3 + $r;
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $x = R(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } return $x; }
my %val;   # value at b (sorted tuple not assumed)
sub tuples { my ($r, $max) = @_; my @out = ([]); for (1..$r) { my @n; for my $t (@out) { my $s = 0; $s += $_ for @$t; for my $x (0..$max - $s) { push @n, [@$t, $x] } } @out = @n; } return @out; }
my @all = tuples($r, $N);
for my $k (0..$h) {
  my %v; for my $b (@all) { $v{join(',', @$b)} = Comps::comp_at($C, $h, $r, $k, [map { $_ + 1 } @$b]); }
  my %alpha;
  for my $a (@all) { my $s = R(0);
    for my $b (@all) { my $ok = 1; for my $i (0..$r-1) { $ok = 0 if $b->[$i] > $a->[$i] } next unless $ok;
      my $w = R(1); my $d = 0; for my $i (0..$r-1) { $w *= binom($a->[$i], $b->[$i]); $d += $a->[$i] - $b->[$i]; }
      $s += (($d % 2) ? -$w : $w) * $v{join(',', @$b)}; }
    $alpha{join(',', @$a)} = $s; }
  print "== genus $h, r = $r, component k = $k ==\n";
  my $low = 0; for my $a (@all) { my $s = 0; $s += $_ for @$a; $low++ if $s <= 2*$h - 3 && $alpha{join(',', @$a)} != 0; }
  print "   non-zero alpha with |a| <= 2h-3: $low\n";
  for my $L (2*$h-2 .. $N) { my @row; for my $a (sort { $a->[0] <=> $b->[0] } grep { my $s = 0; $s += $_ for @$_; $s == $L } @all) { my $x = $alpha{join(',', @$a)}; push @row, "(" . join(',', @$a) . "):$x" if $x != 0; }
    print "   |a| = $L: @row\n"; }
}
