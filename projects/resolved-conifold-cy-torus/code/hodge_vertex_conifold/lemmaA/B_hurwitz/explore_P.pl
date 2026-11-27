#!/usr/bin/perl
# Exploration: the k = 0 polynomial P(eps) = Q(eps - 1) = <Lambda^vee(1) K_{eps-1} lam_R> = sum_n p_n eps^n,
#   p_n = (1/n!) <c(T(-log)) lam_R (psi^2 lam_1)^n>_{h, r+n},  for r = 1 as polynomials in j.
# Usage: perl explore_P.pl h
use strict; use warnings; use lib '.'; use HB; use HBconn; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my $h = shift // 2; my $r = 1; my $n0 = 3*$h - 3 + $r; my $deg = 3*$h - 2;
my %P; # P{j}[n]
for my $j (1..$deg + 2) { my $A = HBconn::A_list([$j], $h, $n0); my $q = HBconn::Q_from_A($A, $h, $r, $j);
  # P(eps) = Q(eps - 1): p_n = sum_m q_m C(m,n) (-1)^(m-n)
  my @p; for my $nn (0..$n0) { my $v = BR(0); for my $m ($nn..$n0) { my $t = $q->[$m] * BR(HB::binom($m, $nn)); $t = -$t if ($m - $nn) % 2; $v += $t } push @p, $v }
  $P{$j} = \@p }
# fit each p_n(j) as a polynomial in j of degree <= deg via Newton interpolation on j = 1..deg+1, check at deg+2
sub fit { my ($xs, $ys) = @_; my $n = @$xs; my @dd = map { $_->copy } @$ys;
  for my $k (1..$n-1) { for (my $i = $n-1; $i >= $k; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($xs->[$i] - $xs->[$i-$k]) } }
  my @acc = (BR(0)) x $n;
  for (my $i = $n-1; $i >= 0; $i--) { my @new = (BR(0)) x $n; for my $k (0..$n-2) { $new[$k+1] += $acc[$k] } for my $k (0..$n-1) { $new[$k] -= $acc[$k] * $xs->[$i] } $new[0] += $dd[$i]; @acc = @new }
  return \@acc }
sub peval { my ($c, $x) = @_; my $v = BR(0); for (my $k = $#$c; $k >= 0; $k--) { $v = $v * $x + $c->[$k] } return $v }
for my $nn (0..$n0) { my @xs = map { BR($_) } 1..$deg + 1; my @ys = map { $P{$_}[$nn] } 1..$deg + 1;
  my $c = fit(\@xs, \@ys); die "fit check failed" unless peval($c, BR($deg + 2)) == $P{$deg + 2}[$nn];
  while (@$c > 1 && $c->[-1] == 0) { pop @$c }
  # integer roots
  my @roots = grep { peval($c, BR($_)) == 0 } -3*$h..4*$h;
  print "p_$nn(j) = ", join(" + ", map { "($c->[$_]) j^$_" } grep { $c->[$_] != 0 } 0..$#$c), "   integer roots: @roots\n" }
