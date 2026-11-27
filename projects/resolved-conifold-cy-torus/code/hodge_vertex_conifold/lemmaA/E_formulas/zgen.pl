# zgen.pl r jmax : all-genus generating polynomials  Zhat_R(u) = Sum_{h>=1} sigma^h Xi_h(R) u^{N_h}/N_h!,  N_h = 2h-3+r,
# for multisets R of integers (parts 0..jmax), using genus 1..6 point evaluations (HodgePt). By Lemma A (proved h <= 6 by the
# engine) only 2h-2 <= D(R) contribute for parts >= 2, so Zhat_R is exact when D(R) <= 10.  Prints the u-coefficients as sigma-polys.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt;
sub R { Math::BigRat->new(shift) }
my ($r, $jmax, $hmax) = @ARGV; $hmax //= 4;
my @pts = ([]); for (1..$r) { my @n; for my $t (@pts) { my $lo = @$t ? $t->[-1] : 0; for my $x ($lo..$jmax) { push @n, [@$t, $x] } } @pts = @n; }
my %val;
for my $h (1..$hmax) { Hodge::setup($h); %HodgePt::memo = ();
  for my $p (@pts) { my $v = HodgePt::WP([], $p); my @s = map { R(0) } 0..$h; for my $k (0..$#$v) { $s[$h - $k] += (($k % 2) ? -1 : 1) * $v->[$k]; } $val{join(',', @$p)}{$h} = \@s; } }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
for my $p (@pts) { my $D = 0; $D += $_ - 1 for @$p; print "R = (", join(',', @$p), "), D = $D:\n";
  for my $h (1..$hmax) { my $N = 2*$h - 3 + $r; my $s = $val{join(',', @$p)}{$h}; my @t; for my $i (0..$#$s) { push @t, "(" . ($s->[$i] / fact($N)) . ")s^$i" if $s->[$i] != 0; }
    print "   u^$N: ", (@t ? join(' + ', @t) : '0'), "\n"; } }
