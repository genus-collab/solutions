#!/bin/bash
# genus0_check.sh : Theorem 1 (xy_npoint.pl, g = 0, P2-normalised) against the genus-0 formula Xi_0(R) = Gamma(J)/Gamma(J-r+3) = (J-1)...(J-r+3).
cd "$(dirname "$0")"
for args in "1 0 3 5" "2 0 3 4" "1 0 4 4" "2 0 4 3" "1 0 5 3" "3 0 5 2"; do
  timeout 900 perl xy_npoint.pl $args /dev/null 2>&1 | perl -ne 'BEGIN{use Math::BigRat; $nt=0; $nb=0;} if (/^tau=(\S+) g=0 mu=\(([\d,]+)\): xy (\S+)/) { my @m = split /,/, $2; my $J = 0; $J += $_ for @m; my $r = @m; my $p = Math::BigRat->new(1); $p *= ($J - $_) for 1..$r-3; $nt++; $nb++ if Math::BigRat->new($3) != $p; $tau=$1; $rr=$r; } END { print "genus 0, tau=$tau, r=$rr: compared $nt with Gamma(J)/Gamma(J-r+3), mismatches $nb\n"; }'
done
