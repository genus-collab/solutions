#!/usr/bin/perl
# Referee (outlook remark): principal-specialisation branching rule used in the referee's comment on the CY obstruction:
#   sum_{mu = nu - box} s_mu(1,q,q^2,...) = (1-q) * sum_{box in nu} q^{c(box)} * s_nu(1,q,q^2,...),
# with s_nu(1,q,q^2,...) = q^{n(nu)} / prod_{box} (1 - q^{h(box)})  (Stanley EC2 7.21.3), checked exactly at q = 1/3, 2/7, |nu| <= 9.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefHur; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
sub conj { my $l = shift; my @c; for my $i (0..$#$l) { $c[$_]++ for 0..$l->[$i] - 1 } \@c }
sub sq { my ($l, $q) = @_; my $c = conj($l); my $n = 0; $n += $_ * $_ for 0..$#$l; my $nn = 0; $nn += $_ * $l->[$_] for 0..$#$l;
  my $v = $q ** $nn; for my $i (0..$#$l) { for my $j (0..$l->[$i] - 1) { my $h = $l->[$i] - $j + $c->[$j] - $i - 1; $v /= (BR(1) - $q ** $h) } } $v }
my ($n, $bad) = (0, 0);
for my $q (BR("1/3"), BR("2/7")) { for my $d (1..9) { for my $nu (@{RefHur::parts_of($d)}) {
  my $lhs = BR(0); for my $i (0..$#$nu) { next if $i < $#$nu && $nu->[$i + 1] == $nu->[$i]; my @m = @$nu; $m[$i]--; @m = grep { $_ > 0 } @m; $lhs += sq(\@m, $q) }
  my $cs = BR(0); for my $i (0..$#$nu) { for my $j (0..$nu->[$i] - 1) { $cs += $q ** ($j - $i) } }
  my $rhs = (BR(1) - $q) * $cs * sq($nu, $q); $n++; $bad++ if $lhs != $rhs } } }
print "branching rule (principal specialisation): $n checks, $bad failures\n";
