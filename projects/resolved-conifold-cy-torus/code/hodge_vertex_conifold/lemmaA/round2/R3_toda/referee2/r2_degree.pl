#!/usr/bin/perl
# Second referee: Route A Prop. 2.1 / Cor. 2.2 (the input I4 of Lemma 2.1) at high genus, from the MV Schur side alone:
#   Q_g(s;B) := sum_m s^m sum_{a<=m} (M_a/a!)(-1)^(m-a) L_{m-a}^{(2g-3+r+a)}(-(J+a)),  M_a = M_g(B+1^a),
# must be a polynomial of degree <= N = 3g-3+r (so that e^{-uJ} G^(g)_B = w^{2g-2+r} Q_g(w-1;B), w = 1+s), and Q_g(-1;B) = Xi_g(B).
# Checks [s^(N+1)] Q = [s^(N+2)] Q = 0 for g < gmax (the unit-leg range of a genus-gmax run covers a <= N_g + 3).
# Usage: perl r2_degree.pl prime tau gmax "B1;B2;..."
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use R2mv;
$| = 1;
my ($prime, $tau, $gmax, $bl) = @ARGV; R2mv::setp($prime); my $p = $prime;
my @Bs = map { [ grep { length } split /,/, $_ ] } split /;/, $bl;
my $res = R2mv::compute($gmax, \@Bs, $tau, report => 1, extra => 3, empty => 1);
my ($n, $bad) = (0, 0);
for my $B (@Bs) { my @b = sort { $b <=> $a } @$B; my $key = join(',', @b); my $r = @b; my $J = 0; $J += $_ for @b;
  for my $g (1..$gmax - 1) { next if $g == 1 && $r == 0; my $N = 3*$g - 3 + $r; my $M = $res->{$key}{M}{$g}; next unless $M && @$M >= $N + 3;
    my @c; for my $m ($N + 1, $N + 2) { my $s = 0; for my $a (0..$m) { next if $r == 0 && $a == 0; my $Ma = $M->[$a]; next unless defined $Ma;
        my $t = R2mv::lag($m - $a, 2*$g - 3 + $r + $a, R2mv::md(-$J - $a)) * $Ma % $p * R2mv::ifac($a) % $p; $s = (($m - $a) & 1) ? R2mv::md($s - $t) : ($s + $t) % $p; } push @c, $s; }
    $n++; my $ok = !($c[0] || $c[1]); $bad++ unless $ok;
    printf "g=%d B=(%s) N=%d: [s^(N+1)]Q = %d, [s^(N+2)]Q = %d  %s\n", $g, join(',', reverse @b), $N, $c[0], $c[1], ($ok ? "OK (deg <= N)" : "FAIL"); } }
printf "SUMMARY degree bound: %d cases, %d failures\n", $n, $bad;
