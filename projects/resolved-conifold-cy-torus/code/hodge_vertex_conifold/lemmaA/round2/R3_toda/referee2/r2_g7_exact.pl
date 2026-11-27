#!/usr/bin/perl
# Second referee: EXACT genus-7 test of Theorem M (P_2) and of the P_3 / P_4 leading-order relations, using only Hodge-integral data:
# genus-7 values from Route E's exact Hodge engine (R3_toda/out/g7_points_A.out: v-polynomials, v = -1/sigma), genus <= 6 values
# from Route E's exact tables, genus 0 closed form.  No Marino-Vafa / Schur input at all.  Each relation divided by tau is a
# polynomial in tau of degree <= 2G + eps + k - 1, so vanishing at 2G + 6 distinct rational points proves the identity in tau.
# Usage: perl r2_g7_exact.pl
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use R2rel; use R2tab;
$| = 1;
my %g7; my $f = "$FindBin::Bin/../out/g7_points_A.out"; open my $fh, '<', $f or die "open $f";
while (<$fh>) { s/[\r\n]//g; next if /^#/; my @x = split /\s*\|\s*/; next unless @x == 3; my @T = sort { $a <=> $b } split /,/, $x[1];
  $g7{join(',', @T)} = [ map { Math::BigRat->new($_) } split /\s*;\s*/, $x[2] ]; } close $fh;
print "# exact genus-7 points: ", join(' ', map { "($_)" } sort keys %g7), "\n";
my $Q = R2rel::ring_Q();
my @taus = map { Math::BigRat->new($_) } qw(1/2 1/3 2/3 3/2 2 3 5/2 -1/3 -2/3 -3/2 -5/2 -3 1/5 4/5 7/3 -7/4 5/7 -9/7 11/5 -13/6 2/9 -11/9);
my $G = 7;
for my $test ([2, [13]], [2, [14]]) { my ($k, $Rp) = @$test; my $P = R2rel::hirota_Pk($k);
  for my $eps (0, 1) { my ($nz, $npts, $missing) = (0, 0, '');
    for my $tau (@taus) { my $sig = $tau * ($tau + 1); my $v = Math::BigRat->new(-1) / $sig;
      my $X = sub { my ($g, $U) = @_; my $ka = join(',', @$U);
        if ($g == 7) { die "MISSING g7 ($ka)\n" unless exists $g7{$ka}; my $c = $g7{$ka}; my $s = Math::BigRat->new(0); my $q = Math::BigRat->new(1); for my $i (0..$#$c) { $s += $c->[$i] * $q; $q *= $v; } return $s; }
        return R2tab::xi_at($g, $U, $tau); };
      my ($val) = eval { R2rel::relation(P => $P, G => $G, eps => $eps, Rparts => $Rp, ring => $Q, tau => $tau, X => $X) };
      if ($@) { $missing = $@; $missing =~ s/\n//g; last; }
      $npts++; $nz++ if $val != 0; }
    my $J = 0; $J += $_ for @$Rp; my $DK = $J - @$Rp + $k - 1;
    printf "P_%d G=%d eps=%d R=(%s) [R+%d: D=%d, excess %d]: %s\n", $k, $G, $eps, join(',', @$Rp), $k, $DK, 2*$G - 2 - $DK,
      ($missing ? "skipped ($missing)" : ($nz ? "NONZERO at $nz of $npts points" : "zero at all $npts rational points (identity in tau)")); } }
