#!/usr/bin/perl
# Referee: independent test of the refuted ansatz C3(f)
#   Xi^(k)_h(R) =? Sum_{g=0}^{h-k} c_g * Xi^(0)_g(R u 0^{2(h-g)}),   c_g independent of R,
# by EVALUATION at lattice points (not by comparing polynomial coefficients as fit_decomp.pl does):
#   target Xi^(k)_h(R) from the referee's own component extraction (out/ucomp_g<h>_r<r>.txt);
#   basis: g = 0: (J-1)...(J-N), N = 2h-3+r (= Xi^(0)_0(R u 0^{2h}));  g >= 1: prod_{i<2(h-g)} (D - i - 2g + 2) * Xi^(0)_g(R),
#   with Xi^(0)_g(R) computed from Hurwitz numbers (ref_elsv.pl), not from the author's tables.
# The linear system over all lattice points with parts >= 1 and J <= Jmax is solved exactly; report consistent / inconsistent.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use rpoly;
require "$FindBin::Bin/ref_elsv.pl";
my %cache; sub X0 { my ($g, @R) = @_; my $k = "$g|@R"; return $cache{$k} //= (main::xi0($g, @R))[0]; }
for my $case ([2, 2, 1, 6], [3, 2, 1, 7], [3, 2, 2, 7], [2, 1, 1, 8], [3, 1, 1, 10]) { my ($h, $r, $k, $Jmax) = @$case;
  open(my $fh, '<', "$FindBin::Bin/out/ucomp_g${h}_r${r}.txt") or die; my $T; while (<$fh>) { s/\r//; chomp; my @f = split /\s*\|\s*/; $T = parse_rp($f[3]) if $f[2] == $k; } close $fh;
  my @pts; if ($r == 1) { @pts = map { [$_] } 1 .. $Jmax; } else { for my $a (1 .. $Jmax) { for my $b ($a .. $Jmax - $a) { push @pts, [$a, $b]; } } }
  my @rows; my $N = 2 * $h - 3 + $r;
  for my $p (@pts) { my $J = 0; $J += $_ for @$p; my $D = $J - $r; my @row;
    my $b0 = Math::BigRat->new(1); $b0 *= ($J - $_) for 1 .. $N; push @row, $b0;
    for my $g (1 .. $h - $k) { my $f = Math::BigRat->new(1); $f *= ($D - $_ - 2 * $g + 2) for 0 .. 2 * ($h - $g) - 1; push @row, $f * X0($g, @$p); }
    push @row, rp_eval($T, @$p); push @rows, \@row; }
  # exact Gaussian elimination on the augmented matrix
  my $nc = $h - $k + 1; my $row = 0; my @M = map { [ map { $_->copy } @$_ ] } @rows;
  for my $col (0 .. $nc - 1) { my $pv = $row; $pv++ while ($pv < @M && $M[$pv][$col]->is_zero); next if $pv >= @M; @M[$row, $pv] = @M[$pv, $row];
    my $inv = Math::BigRat->new(1) / $M[$row][$col]; $M[$row] = [ map { $_ * $inv } @{$M[$row]} ];
    for my $q (0 .. $#M) { next if $q == $row; my $fct = $M[$q][$col]; next if $fct->is_zero; $M[$q] = [ map { $M[$q][$_] - $fct * $M[$row][$_] } 0 .. $nc ]; } $row++; }
  my $cons = 1; for my $q ($row .. $#M) { $cons = 0 unless $M[$q][$nc]->is_zero; }
  printf "genus %d r %d k %d: %d lattice points, %d basis functions: %s\n", $h, $r, $k, scalar(@pts), $nc,
    ($cons ? "CONSISTENT (c_g = " . join(', ', map { $M[$_][$nc]->bstr } 0 .. $row - 1) . ")" : "INCONSISTENT -> ansatz refuted"); }
