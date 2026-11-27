#!/usr/bin/perl
# Referee check of the literature input (I2)-(I3) for GENERIC content-product 2D Toda tau functions at fixed second times:
#   tau_n(p) = Sum_{|nu| <= d} B^{f_2(nu) + n|nu|} s_nu(p') s_nu(p),   B, p'_1..p'_d random rationals,
# (MV: B = e^{tau x}, p' = q^rho, so tau_0 = Z_MV, tau_1 = Z_MV(e^{tau x} o p)).  Exact (Math::BigRat), all monomials of degree <= d - k.
# Tested Hirota equations (t_k = p_k/k):
#   (JM1)  (D_1^2 + D_2) tau_0.tau_1 = 0              [= 2 P_2; the author's first mKP equation]
#   (JM2)  (D_1^3 + 3 D_1 D_2 - 4 D_3) tau_1.tau_0 = 0  [Jimbo-Miwa's second mKP equation]
#   (Pk)   (D_k + p_k(-D~)) tau_0.tau_1 = 0, k = 2..kmax  [the author's u_k-linear members]
# Controls (must fail): the same with tau_1 replaced by tau_{-1}, and with tau_0, tau_1 swapped.
# Usage: perl ref_toda_generic.pl d kmax B seed
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use ref_lib;
$| = 1;
sub R { Math::BigRat->new(@_) }
my ($d, $kmax, $Bs, $seed) = @ARGV; $d //= 7; $kmax //= 4; $Bs //= "3/2"; $seed //= 1;
my $B = R($Bs); srand($seed);
my @pp = (undef, map { R(int(rand(19)) - 9) / R(1 + int(rand(7))) } 1..$d);   # p'_k
print "B = $B, p' = (", join(', ', map { "$_" } @pp[1..$d]), ")\n";
# Schur functions
my %S; my %Sval;
for my $n (0..$d) { for my $nu (ref_lib::partitions($n)) { my $kn = join(',', @$nu); my %s; my $val = R(0);
    for my $mu (ref_lib::partitions($n)) { my $c = ref_lib::chi($nu, $mu); next unless $c; my $z = ref_lib::zee(@$mu);
      $s{ref_lib::mkey(@$mu)} = R($c) / $z; my $pv = R(1); $pv *= $pp[$_] for @$mu; $val += R($c) * $pv / $z; }
    $S{$kn} = \%s; $Sval{$kn} = $val; } }
# sanity: s_(n) = h_n, s_(1^n) = e_n via Newton (check at n = 3): s_(2,1) = (p1^3 - p3)/3
{ my $s21 = $S{"2,1"}; my $ok = ($s21->{"1,1,1"} == R(1,3)) && ($s21->{"3"} == R(-1,3)) && !exists $s21->{"2,1"};
  print "sanity s_(2,1) = (p1^3 - p3)/3: ", ($ok ? "OK" : "FAIL"), "\n"; }
sub tau { my ($n) = @_; my $t = {};
  for my $kn (keys %S) { my @nu = $kn eq '' ? () : split /,/, $kn; my $sz = 0; $sz += $_ for @nu;
    my $w = $B ** (ref_lib::f2(@nu) + $n * $sz) * $Sval{$kn}; next if $w == 0; $t = ref_lib::padd($t, $S{$kn}, $w); }
  return $t; }
my %T = map { $_ => tau($_) } (-1, 0, 1);
# Hirota polynomial: list of [coef, [indices]]
sub apply { my ($P, $f, $g, $wt) = @_; my $res = {};
  for my $t (@$P) { my ($c, $idx) = @$t; $res = ref_lib::padd($res, ref_lib::hirota($f, $g, $idx, $d - $wt), $c); } return $res; }
sub report { my ($name, $res) = @_; my @nz = grep { $res->{$_} != 0 } keys %$res;
  printf "%-46s %s\n", $name, (@nz ? "NONZERO in " . scalar(@nz) . " monomials (e.g. p_[" . $nz[0] . "]: " . $res->{$nz[0]} . ")" : "0 (all monomials of degree <= " . ($d) . " - weight)"); return scalar @nz; }
my $JM1 = [[R(1), [1,1]], [R(1), [2]]];
my $JM2 = [[R(1), [1,1,1]], [R(3), [1,2]], [R(-4), [3]]];
report("(JM1) (D1^2+D2) tau0.tau1", apply($JM1, $T{0}, $T{1}, 2));
report("  control (D1^2+D2) tau1.tau0", apply($JM1, $T{1}, $T{0}, 2));
report("  control (D1^2+D2) tau0.tau_{-1}", apply($JM1, $T{0}, $T{-1}, 2));
report("  control (D1^2+D2) tau_{-1}.tau0", apply($JM1, $T{-1}, $T{0}, 2));
report("(JM2) (D1^3+3D1D2-4D3) tau1.tau0", apply($JM2, $T{1}, $T{0}, 3));
report("  control (JM2) tau0.tau1", apply($JM2, $T{0}, $T{1}, 3));
for my $k (2..$kmax) { my @P = ([R(1), [$k]], ref_lib::pk_minus_Dtilde($k));
  report("(P$k) (D_$k + p_$k(-D~)) tau0.tau1", apply(\@P, $T{0}, $T{1}, $k));
  report("  control (P$k) tau1.tau0", apply(\@P, $T{1}, $T{0}, $k));
  report("  control (P$k) tau0.tau_{-1}", apply(\@P, $T{0}, $T{-1}, $k)); }
# KP for a single tau (sanity of the Hirota code): (D1^4 + 3 D2^2 - 4 D1 D3) tau.tau = 0
report("(KP) (D1^4+3D2^2-4D1D3) tau0.tau0", apply([[R(1),[1,1,1,1]], [R(3),[2,2]], [R(-4),[1,3]]], $T{0}, $T{0}, 4));
report("(KP) (D1^4+3D2^2-4D1D3) tau1.tau1", apply([[R(1),[1,1,1,1]], [R(3),[2,2]], [R(-4),[1,3]]], $T{1}, $T{1}, 4));
