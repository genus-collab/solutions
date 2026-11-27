#!/usr/bin/perl
# (B2) The finite formula for the Hurwitz component, checked against the exact Hodge data.
#   Xi^(0)_h(R) = sum_{c=0}^{n0} (-1)^c (A_h(R,1^c)/c!) sum_{i=0}^{n0-c} C(n0+A0, n0-c-i) (J+c)^i / i!,
#   A0 = 2h-2+r, n0 = 3h-3+r, A_h(R,1^c) from ELSV (connected Hurwitz numbers H_{h, R+1^c}).
# Also the whole s-polynomial Q(s) = <Lambda^vee(1) K_s lam_R> (K_s = exp(sum((1+s)^k-1)kappa_k/k)) against q_g<g>_r<r>.txt.
# Usage: perl b2_check.pl
use strict; use warnings; use lib '.'; use HB; use HBconn; use Data; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my ($n, $bad) = (0, 0);
sub report { my ($what, $a, $b) = @_; $n++; my $ok = ($a == $b) ? "ok" : "MISMATCH"; $bad++ if $ok ne "ok"; print "$what: Hurwitz $a, data $b  $ok\n" }
# 1. q tables: full s-polynomial
for my $spec (['q_g2_r1.txt', 2, 1], ['q_g2_r2.txt', 2, 2], ['q_g3_r1.txt', 3, 1], ['q_g3_r2.txt', 3, 2]) {
  my ($file, $h, $r) = @$spec; next unless -e "../data/$file";
  my $ast = Data::load($file)->[0]{ast}; my $n0 = 3*$h - 3 + $r;
  my @tuples = $r == 1 ? (map { [$_] } 1..($h == 2 ? 6 : 8)) : (map { my $a = $_; map { [$a, $_] } $a..($h == 2 ? 6 : 5) } 1..($h == 2 ? 6 : 5));
  for my $R (@tuples) {
    my $J = 0; $J += $_ for @$R;
    # data: Hurwitz component of Q(s) at s = 0..n0 (polynomial of degree <= n0 in s)
    my @qd;
    for my $sv (0..$n0) { my %env = (s => BR($sv)); $env{"j$_"} = BR($R->[$_ - 1]) for 1..$r;
      push @qd, Data::components($ast, $h, $r, \%env)->[0] }
    my @Rn = grep { $_ >= 2 } @$R; my $c1 = grep { $_ == 1 } @$R;
    # Hurwitz side; parts equal to 1 among the probes are treated as unit parts (A_h(R,1^c) with the probe units included)
    my $A = HBconn::A_list(\@Rn, $h, $n0 + $c1);
    # A for the probe multiset R (with c1 probe units) and c extra units: A_h(Rn, 1^(c1+c)); combinatorics: the probes are labelled,
    # the extra units are not; the formula Q_from_A needs A_h(R,1^c) for the probe set R, i.e. A_h(Rn,1^(c1+c))
    my @AR = map { $A->[$c1 + $_] } 0..$n0;
    my $q = HBconn::Q_from_A(\@AR, $h, $r, $J);
    for my $sv (0..$n0) { my $val = BR(0); for my $m (0..$n0) { $val += $q->[$m] * BR($sv) ** $m }
      report("Q h=$h R=(@$R) s=$sv", $val, $qd[$sv]) }
    my $xi = HBconn::xi0_from_A(\@AR, $h, $r, $J);
    my %env = (s => BR(-1)); $env{"j$_"} = BR($R->[$_ - 1]) for 1..$r;
    report("Xi0 h=$h R=(@$R) (s=-1)", $xi, Data::components($ast, $h, $r, \%env)->[0]) } }
# 2. r = 1 closed polynomials from the stored Hodge data (hvc_which_class_out.txt, check_closed_formula_r1.pl), h = 2, 3, 4
sub pr { my ($j, $h) = @_; my $p = BR(1); $p *= BR($j - $_) for 1..(2*$h - 2); return $p }
my %closed = (
  2 => sub { my $j = BR(shift); pr($j, 2) * ($j - 4) * (5*$j - 6) / 5760 },
  3 => sub { my $j = BR(shift); pr($j, 3) * ($j - 6) * (35*$j*$j - 126*$j + 80) / 2903040 },
  4 => sub { my $j = BR(shift); pr($j, 4) * (8064 - 19856*$j + 12436*$j**2 - 2660*$j**3 + 175*$j**4) / 1393459200 } );
for my $h (2..4) { for my $j (1..(3*$h - 1)) { report("Xi0_$h($j) vs stored r=1 polynomial", HBconn::xi0([$j], $h), $closed{$h}->($j)) } }
print "B2: $n comparisons, $bad mismatches\n";
