#!/usr/bin/perl
# t_kpmod.pl -- cross-check of the modular engine kpmod.pm against the exact rows of kp_point.pl (out/beta_tau2_7_K10.txt)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod; use Math::BigInt;
my $tau = mfrac(2, 7); my $S = kp_setup($tau, 1, 10); my $rows = kp_rows($S);
my %ex; open my $fh, '<', "$FindBin::Bin/out/beta_tau2_7_K10.txt" or die; while (<$fh>) { if (/^REL (\d+): (.*)$/) { $ex{$1} = [ split ' ', $2 ]; } } close $fh;
my ($n, $bad) = (0, 0);
for my $r (@$rows) { my ($k, $eq) = @$r; for my $i (0..$#$eq) { my $q = $ex{$k}[$i]; my ($a, $b) = split m{/}, $q; $b //= 1;
    my $am = Math::BigInt->new($a)->bmod($P)->numify; my $bm = Math::BigInt->new($b)->bmod($P)->numify; my $want = mm($am, mi($bm));
    $n++; $bad++ if $want != $eq->[$i]; } }
print "t_kpmod: $n entries compared, $bad mismatches\n";
