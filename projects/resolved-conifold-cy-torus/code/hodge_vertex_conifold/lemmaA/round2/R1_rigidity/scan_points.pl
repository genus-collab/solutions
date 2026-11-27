#!/usr/bin/perl
# scan_points.pl -- pointwise forcing at the special framings tau0 = -k/a (gcd(k,a) = 1, k <= a/2), modular.
# For every genus G <= Gmax and every odd level D in [3, 2G-3] (single-level problem of (OL): y_L = 0 for L < D),
# at tau0 the values y_L(tau0) = Pi_{L+1} Xi_G(L+1) (L = D..2G-3) satisfy
#   (i)   y_L(tau0) = 0 if Pi_{L+1}(tau0) = 0                                   [one-leaf regularity]
#   (ii)  Y(T)(tau0) = 0 if Pi_T(tau0) = 0, |T| = 2 (mode two) or all T (full)  [multi-leaf regularity, first order]
#   (iii) the closed-form relations of orders D..2G-3 (cf = 1)
# and we test whether y_D(tau0) is forced to vanish (rank test: column 0 not in the span of the other columns).
# If a does not divide D+1, forcing gives P_{D+1}(sigma0) = 0 for P_{D+1} = sigma^G Xi_G(D+1), sigma0 = -k(a-k)/a^2.
# Usage: perl scan_points.pl Gmax mode cf [amax]      output: one line per (tau0, G, D) plus a summary per (G, D)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod;
my ($Gmax, $mode, $cf, $amax) = @ARGV; $mode //= 'two'; $cf //= 1; my $Lmax = 2*$Gmax - 3; $amax //= 2*$Gmax - 2;
sub gcd { my ($a, $b) = @_; ($a, $b) = ($b, $a % $b) while $b; return $a; }
sub rank_cols { my ($rows, $cols) = @_; my @M = map { my $r = $_; [ map { $r->[$_] } @$cols ] } @$rows; my $nc = @$cols; my $rk = 0;
  for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#M) { if ($M[$i][$c]) { $pr = $i; last; } } next unless defined $pr; @M[$rk, $pr] = @M[$pr, $rk];
    my $inv = mi($M[$rk][$c]); $M[$rk] = [ map { mm($_, $inv) } @{$M[$rk]} ];
    for my $i (0..$#M) { next if $i == $rk || !$M[$i][$c]; my $f = $M[$i][$c]; for my $j ($c..$nc-1) { $M[$i][$j] = md($M[$i][$j] - mm($f, $M[$rk][$j])) if $M[$rk][$j]; } }
    $rk++; } return $rk; }
my %count;   # count{G,D} = number of sigma0 zeros found
my %pts;
$| = 1;
for my $a (2..$amax) { for my $k (1..int($a/2)) { next unless gcd($k, $a) == 1;
    my $tau = mfrac(-$k, $a); my $S = kp_setup($tau, 1, $Lmax); my $rows = $cf ? kp_rows($S) : [];
    my $Pi = $S->{Pi};
    for my $G (2..$Gmax) { my $Lg = 2*$G - 3;
      for (my $D = 3; $D <= $Lg; $D += 2) {
        my @cols = map { $_ - 1 } ($D..$Lg);   # column index in S (D=1 based): L -> L-1
        my @eqs;
        for my $L ($D..$Lg) { push @eqs, $S->{Y}{$L+1} if !$Pi->($L+1); }
        if ($mode ne 'none') { for my $T (@{$S->{allT}}) { next if $mode eq 'two' && @$T != 2; my $lv = 0; $lv += $_ - 1 for @$T; next if $lv < $D || $lv > $Lg;
            my $p = 1; $p = mm($p, $Pi->($_)) for @$T; push @eqs, $S->{Y}{join(',', @$T)} unless $p; } }
        if ($cf) { for my $r (@$rows) { my ($ko, $eq) = @$r; next if $ko < $D || $ko > $Lg; push @eqs, $eq; } }
        # restrict to levels >= D (columns), and only rows with support there (lower y's are zero by hypothesis)
        my $r1 = rank_cols(\@eqs, \@cols); my @c2 = @cols[1..$#cols]; my $r2 = @c2 ? rank_cols(\@eqs, \@c2) : 0;
        my $forced = ($r1 - $r2 == 1) ? 1 : 0; my $useful = (($D+1) % $a) ? 1 : 0;
        $count{"$G,$D"} += $forced && $useful; push @{$pts{"$G,$D"}}, "$k/$a" if $forced && $useful;
        printf "POINT a=%d k=%d G=%d D=%d forced=%d useful=%d rank=%d/%d\n", $a, $k, $G, $D, $forced, $useful, $r1, scalar(@cols);
      } } } }
for my $G (2..$Gmax) { for (my $D = 3; $D <= 2*$G-3; $D += 2) { my $c = $count{"$G,$D"} // 0;
    printf "SUMMARY G=%d D=%d zeros=%d needed=%d %s  points: %s\n", $G, $D, $c, $G-1, ($c >= $G-1 ? "CLOSES" : "short"), join(' ', @{$pts{"$G,$D"} // []}); } }
