package YFun;
# The pull-back functional  Y_h(j_1..j_r) := Int_{Mbar_{h,r}} pi^*(Omega K_h) prod_i psi_i/(1 - j_i psi_i)
#   = << Sum_{sigma in S_r} prod_{cycles c} kappa(j_c) >>,   kappa(j_c) = Sum_b h_b(j_c) kappa_b,   << X >> = Int_{Mbar_h} Omega K_h X,
# evaluated from the engine's T-table (Hodge::T), as sigma-polynomials of sigma^h Y (sigma^h Y = Sum_k (-1)^k sigma^{h-k} Y^(k)).
use strict; use warnings; use Math::BigRat; use Hodge;
sub R { Math::BigRat->new(shift) }
# set partitions of 0..r-1
sub setparts { my $r = shift; return ([]) if $r == 0; my @out;
  for my $p (setparts($r - 1)) { for my $i (0..$#$p) { my @q = map { [@$_] } @$p; push @{$q[$i]}, $r - 1; push @out, \@q; } push @out, [(map { [@$_] } @$p), [$r - 1]]; }
  return @out; }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
my %Ymemo;
# Ypsi(a) = Int pi^*(Omega K_h) prod psi_i^{a_i}  (a_i >= 1), as v-poly (v = -1/sigma)
sub Ypsi { my @a = sort { $a <=> $b } @_; my $key = join(',', @a); return $Ymemo{$key} if exists $Ymemo{$key};
  my $h = $Hodge::HCUR; my $res = [R(0)];
  for my $P (setparts(scalar @a)) { my $w = R(1); my @beta; my $scal = R(1);
    for my $B (@$P) { $w *= fact(scalar(@$B) - 1); my $b = 0; $b += $a[$_] - 1 for @$B; if ($b == 0) { $scal *= (2*$h - 2) } else { push @beta, $b } }
    my $t = Hodge::T(@beta); $res = Hodge::vadd($res, Hodge::vscale($t, $w * $scal)); }
  $Ymemo{$key} = $res; return $res; }
# sigma^h Y_h(j) at an integer/rational point j, as sigma-poly (array, index = power of sigma)
sub Ysig { my @j = @_; my $h = $Hodge::HCUR; my $r = scalar @j; my $D = 3*$h - 3; my $tot = [map { R(0) } 0..$h];
  my @tuples = ([]); for (1..$r) { my @n; for my $t (@tuples) { my $s = 0; $s += $_ - 1 for @$t; for my $x (1..($D - $s + 1)) { push @n, [@$t, $x] } } @tuples = @n; }
  for my $a (@tuples) { my $mon = R(1); for my $i (0..$r-1) { $mon *= R($j[$i]) ** ($a->[$i] - 1) } next if $mon == 0;
    my $v = Ypsi(@$a); for my $k (0..$#$v) { $tot->[$h - $k] += (($k % 2) ? -1 : 1) * $v->[$k] * $mon; } }
  return $tot; }
1;
