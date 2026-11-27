# hessian.pl h : the Hessian of proof_M_residue.md s.8,
#   H(j1,j2) = <Omega K Pt_{j1} Pt_{j2}>,  Pt_j = lam_j - j lam_1 + (j-1)  (kills j = 0, 1),
#   H = Xi(j1,j2) - [j1 + (2h-1)(j2-1)] Xi(j1) - [j2 + (2h-1)(j1-1)] Xi(j2) + (2h-1)(2h-2)(j1-1)(j2-1) F_h(q_0),  F_h(q_0) = Xi(0)/(2-2h).
# Checks divisibility by j1(j1-1) j2(j2-1), and whether the quotient vanishes on lines J = l or j_i = l (per component).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use CPoly;
sub R { Math::BigRat->new(shift) }
my $h = shift;
my $C2 = CPoly::load($FindBin::Bin, $h, 2); my $C1 = CPoly::load($FindBin::Bin, $h, 1);
sub divlin { my ($p, $v, $a) = @_; my $d = Poly::degree_in($p, $v); return {} if $d < 0; my @c = map { Poly::coeff_of($p, $v, $_) } 0..$d;
  my @q; my $carry = {}; for (my $k = $d; $k >= 1; $k--) { $carry = Poly::add($c[$k], Poly::scale($carry, $a)); $q[$k-1] = $carry; }
  my $rem = Poly::add($c[0], Poly::scale($carry, $a)); return undef unless Poly::is_zero($rem);
  my $res = {}; for my $k (0..$#q) { $res = Poly::add($res, Poly::mul($q[$k], Poly::pow(Poly::var($v), $k))); } return $res; }
for my $k (0..$h) {
  my $X2 = $C2->[$k]; my $X1a = $C1->[$k]; my $X1b = CPoly::subst_poly($X1a, { j1 => Poly::var('j2') });
  my $F = Poly::scale(CPoly::subst_poly($X1a, { j1 => Poly::const(0) }), R(1) / (2 - 2*$h));
  my $j1 = Poly::var('j1'); my $j2 = Poly::var('j2'); my $one = Poly::const(1);
  my $c1 = Poly::add($j1, Poly::scale(Poly::sub_($j2, $one), 2*$h - 1)); my $c2 = Poly::add($j2, Poly::scale(Poly::sub_($j1, $one), 2*$h - 1));
  my $H = Poly::sub_($X2, Poly::add(Poly::mul($c1, $X1a), Poly::mul($c2, $X1b)));
  $H = Poly::add($H, Poly::mul(Poly::scale(Poly::mul(Poly::sub_($j1, $one), Poly::sub_($j2, $one)), (2*$h-1)*(2*$h-2)), $F));
  my $q = $H; my $ok = 1; for my $v ('j1', 'j2') { for my $a (0, 1) { $q = divlin($q, $v, $a); unless (defined $q) { $ok = 0; last } } last unless $ok; }
  unless ($ok) { print "genus $h k=$k: H not divisible by j1(j1-1)j2(j2-1)?!\n"; next; }
  my @lines; for my $l (-4..3*$h) { my $t = CPoly::subst_poly($q, { j2 => Poly::sub_(Poly::const($l), Poly::var('j1')) }); push @lines, $l if Poly::is_zero($t); }
  my @pl; for my $l (-4..3*$h) { my $t = CPoly::subst_poly($q, { j1 => Poly::const($l) }); push @pl, $l if Poly::is_zero($t); }
  print "genus $h k=$k: H = j1(j1-1)j2(j2-1) Q, deg Q = ", Poly::total_degree($q, 'j1', 'j2'), "; Q vanishes on J = {@lines}; on j1 = {@pl}\n"; }
