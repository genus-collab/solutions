#!/usr/bin/perl
# Second referee: validate the from-scratch MV Schur-side engine (R2mv.pm) against Route E's exact Hodge tables (genus <= 6)
# and against the exact genus-7 Hodge-engine points of R3_toda/out/g7_points_A.out (v-polynomials, v = -1/sigma).
# Usage: perl r2_validate.pl prime tau_mod_p gmax "B1;B2;..."     (B = comma-separated parts >= 2)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use R2mv; use R2tab;
$| = 1;
my ($prime, $taup, $gmax, $blist) = @ARGV; R2mv::setp($prime); my $p = $prime;
my @Bs = map { [split /,/, $_] } split /;/, $blist;
my $t0 = time; my $res = R2mv::compute($gmax, \@Bs, $taup, report => 1, verbose => ($ENV{VERB} // 0));
printf "# engine done in %d s (p = %d, tau = %d)\n", time - $t0, $p, $taup;
my $sig = $taup * R2mv::md($taup + 1) % $p; my $v = R2mv::md(-R2mv::inv($sig));
# genus-7 exact engine points
my %g7; my $g7f = "$FindBin::Bin/../out/g7_points_A.out";
if (-e $g7f) { open my $fh, '<', $g7f or die; while (<$fh>) { s/[\r\n]//g; next if /^#/; my @f = split /\s*\|\s*/; next unless @f == 3; my @T = sort { $a <=> $b } split /,/, $f[1];
    $g7{"$f[0]|" . join(',', @T)} = [ map { Math::BigRat->new($_) } split /\s*;\s*/, $f[2] ]; } close $fh; }
my ($nok, $nbad, $nskip) = (0, 0, 0);
for my $B (@Bs) { my @T = sort { $a <=> $b } @$B; my $key = join(',', sort { $b <=> $a } @$B); my $r = @T; my $J = 0; $J += $_ for @T; my $D = $J - $r;
  for my $g (1..$gmax) { my $mine = $res->{$key}{Xi}{$g}; next unless defined $mine; my $ref; my $src;
    if ($g <= 6 && R2tab::available($g, $r)) { my $c = R2tab::comps($g, \@T); my $s = 0; my $q = 1; for my $k (0..$#$c) { $s = ($s + R2mv::rat($c->[$k]) * $q) % $p; $q = $q * $v % $p; } $ref = $s; $src = "RouteE-table"; }
    elsif (exists $g7{"$g|" . join(',', @T)}) { my $c = $g7{"$g|" . join(',', @T)}; my $s = 0; my $q = 1; for my $k (0..$#$c) { $s = ($s + R2mv::rat($c->[$k]) * $q) % $p; $q = $q * $v % $p; } $ref = $s; $src = "g7-engine"; }
    my $exc = 2*$g - 2 - $D;
    if (defined $ref) { my $ok = ($ref == $mine); $ok ? $nok++ : $nbad++;
      printf "g=%d T=(%s) D=%d excess=%d: R2mv %d, %s %d  %s\n", $g, join(',', @T), $D, $exc, $mine, $src, $ref, ($ok ? "MATCH" : "MISMATCH"); }
    else { $nskip++; printf "g=%d T=(%s) D=%d excess=%d: R2mv %d (no reference)%s\n", $g, join(',', @T), $D, $exc, $mine, ($exc >= 1 ? ($mine == 0 ? "  [in range: 0 as Lemma A predicts]" : "  [in range: NONZERO]") : ""); } } }
printf "SUMMARY validate p=%d tau=%d gmax=%d: %d matches, %d mismatches, %d without reference; %d s\n", $p, $taup, $gmax, $nok, $nbad, $nskip, time - $t0;
