#!/usr/bin/perl
# Second referee: high-genus tests (G = 7..10) of the leading-order relations of the u_k-linear members P_k of the 2D Toda / mKP
# hierarchy (Theorem M for k = 2, Theorem M_k for k = 3, 4), outside the Lemma-A range where the relations are non-trivial,
# with ALL genus >= 1 values computed by the referee's own MV Schur-side engine R2mv.pm (mod p), and the genus <= 6 values
# cross-checked against Route E's exact tables. Also direct Lemma-A checks (values in the range must be 0) and perturbation controls.
# Usage: perl r2_highg.pl prime seed gmax "tests" "inrange"
#   tests   = ";"-separated  k:G:parts   e.g. "2:10:19;2:10:20;3:10:19"   (parts comma-separated, R may be empty: "2:7:")
#   inrange = ";"-separated  profiles to evaluate at genus gmax, e.g. "18;2,17"
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use R2mv; use R2rel; use R2tab;
$| = 1;
my ($prime, $seed, $gmax, $tests, $inr) = @ARGV; $inr //= "";
R2mv::setp($prime); my $p = $prime; srand($seed); my $tau = 2 + int(rand($p - 10));
my $ring = R2rel::ring_P($p);
my @tests = map { my ($k, $G, $pp) = split /:/, $_, 3; [ $k, $G, [ grep { length } split /,/, ($pp // '') ] ] } grep { length } split /;/, $tests;
my @inr = map { [ split /,/, $_ ] } grep { length } split /;/, $inr;
printf "# referee-2 high-genus test: p = %d, tau = %d (seed %d), gmax = %d\n", $p, $tau, $seed, $gmax;
# 1. dry run: which (g, U) does each relation need?
my %need;
for my $t (@tests) { my ($k, $G, $Rp) = @$t; my $P = R2rel::hirota_Pk($k);
  for my $eps (0, 1) { R2rel::relation(P => $P, G => $G, eps => $eps, Rparts => $Rp, ring => $ring, tau => $tau,
      X => sub { my ($g, $U) = @_; $need{join(',', @$U)}{$g} = 1 if $g >= 1; 1 }); } }
for my $U (@inr) { $need{join(',', sort { $a <=> $b } @$U)}{$gmax} = 1; }
my @Ulist = map { [ split /,/, $_ ] } sort keys %need;
printf "# %d profiles needed: %s\n", scalar(@Ulist), join(' ', map { "(" . join(',', @$_) . ")" } @Ulist);
# 2. MV engine
my $t0 = time; my $res = R2mv::compute($gmax, \@Ulist, $tau, report => 1);
printf "# MV engine done in %d s\n", time - $t0;
my $sig = $tau * R2mv::md($tau + 1) % $p; my $vv = R2mv::md(-R2mv::inv($sig));
my %X; my ($ntab, $badtab) = (0, 0);
for my $U (@Ulist) { my $kd = join(',', sort { $b <=> $a } @$U); my $ka = join(',', sort { $a <=> $b } @$U);
  for my $g (1..$gmax) { my $v = $res->{$kd}{Xi}{$g}; next unless defined $v; $X{"$g|$ka"} = $v;
    if ($g <= 6 && R2tab::available($g, scalar @$U)) { my $c = R2tab::comps($g, [sort { $a <=> $b } @$U]); my $s = 0; my $q = 1; for my $k (0..$#$c) { $s = ($s + R2mv::rat($c->[$k]) * $q) % $p; $q = $q * $vv % $p; }
      $ntab++; if ($s != $v) { $badtab++; print "# TABLE MISMATCH g=$g U=($ka): engine $v table $s\n"; } } } }
printf "# cross-check against Route E tables (genus <= 6): %d values, %d mismatches\n", $ntab, $badtab;
my %pert;
my $Xcb = sub { my ($g, $U) = @_; my $ka = join(',', @$U); if ($g == 0) { return R2mv::rat(R2tab::xi0($U)); }
  my $key = "$g|$ka"; die "missing value $key" unless exists $X{$key}; my $v = $X{$key}; $v = ($v + $pert{$key}) % $p if exists $pert{$key}; $v };
# 3. relations
my ($nrel, $nbad) = (0, 0);
for my $t (@tests) { my ($k, $G, $Rp) = @$t; my $P = R2rel::hirota_Pk($k); my $J = 0; $J += $_ for @$Rp; my $DK = $J - @$Rp + $k - 1;
  for my $eps (0, 1) { my ($v, $nt) = R2rel::relation(P => $P, G => $G, eps => $eps, Rparts => $Rp, ring => $ring, tau => $tau, X => $Xcb);
    $nrel++; $nbad++ if $v;
    printf "P_%d G=%d eps=%d R=(%s) [R+%d: D=%d, excess %d] %d terms: residual %d %s\n", $k, $G, $eps, join(',', @$Rp), $k, $DK, 2*$G - 2 - $DK, $nt, $v, ($v ? "FAIL" : "OK"); }
  # control: perturb the genus-G target value X_G(R + k)
  my $tk = join(',', sort { $a <=> $b } (@$Rp, $k)); %pert = ("$G|$tk" => 1);
  my $vc = R2rel::relation(P => $P, G => $G, eps => 0, Rparts => $Rp, ring => $ring, tau => $tau, X => $Xcb); %pert = ();
  printf "   control (X_%d(%s) + 1, eps=0): residual %d %s\n", $G, $tk, $vc, ($vc ? "(nonzero, as it must be)" : "ZERO - control failed"); }
# 4. direct Lemma A values in the range
for my $U (@inr) { my $ka = join(',', sort { $a <=> $b } @$U); my $r = @$U; my $J = 0; $J += $_ for @$U; my $D = $J - $r;
  my $v = $X{"$gmax|$ka"}; printf "Lemma A direct: Xi_%d(%s), D = %d (range D <= %d): %d %s\n", $gmax, $ka, $D, 2*$gmax - 3, $v, ($D <= 2*$gmax - 3 ? ($v ? "NONZERO IN RANGE" : "= 0 OK") : ($v ? "(outside range, nonzero)" : "(outside range, zero)")); }
printf "SUMMARY highg p=%d gmax=%d: %d relations, %d failures; table cross-check %d/%d; total %d s\n", $p, $gmax, $nrel, $nbad, $ntab - $badtab, $ntab, time - $t0;
