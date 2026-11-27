#!/usr/bin/perl
# Test of the decomposition hypothesis (exploratory, route C3):
#   Xi^(k)_h(R) =? Sum_{g=0}^{h-k} c_{h,k,g} * Xi^H_g(R u 0^{2(h-g)})
# with Xi^H_g(R u 0^m) = prod_{i=0}^{m-1} (D_R - i - 2g + 2) * Xi^H_g(R) (dilaton relation, g >= 1) and
# Xi^H_0(R u 0^{2h}) = (J-1)(J-2)...(J-N), N = 2h-3+r (lambda_g formula divided by b_h).
# Reads comp_g<g>_r<r>.txt (or _psi.txt) for the Hurwitz components Xi^H_g = component k = 0.
# Usage: perl fit_decomp.pl h r
use strict; use warnings;
use FindBin; use lib $FindBin::Bin;
use polylib; use Math::BigRat;
my ($h, $r) = @ARGV;
sub readcomp { my ($g, $r) = @_; for my $f ("comp_g${g}_r${r}.txt", "comp_g${g}_r${r}_psi.txt") { my $p = "$FindBin::Bin/$f"; next unless -e $p;
  open(my $fh, '<', $p) or die; my @X; while (<$fh>) { chomp; my @f = split /\s*\|\s*/; $X[$f[2]] = parse_expr($f[3], {}); } close $fh; return \@X; } die "no comp file for g=$g r=$r"; }
my $Xh = readcomp($h, $r);
my $J = {}; $J = padd($J, { "j$_" => Math::BigRat->new(1) }) for 1 .. $r;
my $Dp = padd($J, pconst(-$r));
my @basis;
for my $g (0 .. $h) { my $b;
  if ($g == 0) { my $N = 2 * $h - 3 + $r; $b = pconst(1); for my $l (1 .. $N) { $b = pmul($b, padd($J, pconst(-$l))); } }
  else { my $Xg = readcomp($g, $r)->[0]; $b = $Xg; for my $i (0 .. 2 * ($h - $g) - 1) { $b = pmul($b, padd($Dp, pconst(-$i - 2 * $g + 2))); } }
  push @basis, $b; }
for my $k (1 .. $h - 1) {
  my @B = @basis[0 .. $h - $k]; my $T = $Xh->[$k];
  # linear system: Sum_g c_g B_g = T, coefficientwise in the monomials
  my %keys; $keys{$_} = 1 for (keys %$T, map { keys %$_ } @B); my @ks = sort keys %keys;
  my $nc = scalar @B; my @M; for my $kk (@ks) { push @M, [ (map { pcoef($_, $kk) } @B), pcoef($T, $kk) ]; }
  # Gaussian elimination
  my $row = 0; my @pivcol;
  for my $col (0 .. $nc - 1) { my $p = $row; $p++ while ($p < @M && $M[$p][$col]->is_zero); next if $p >= @M; @M[$row, $p] = @M[$p, $row];
    my $inv = Math::BigRat->new(1) / $M[$row][$col]; $M[$row] = [ map { $_ * $inv } @{$M[$row]} ];
    for my $q (0 .. $#M) { next if $q == $row; my $f = $M[$q][$col]; next if $f->is_zero; $M[$q] = [ map { $M[$q][$_] - $f * $M[$row][$_] } 0 .. $nc ]; }
    push @pivcol, $col; $row++; }
  my $consistent = 1; for my $q ($row .. $#M) { $consistent = 0 unless $M[$q][$nc]->is_zero; }
  my @c = map { Math::BigRat->new(0) } 0 .. $nc - 1; for my $i (0 .. $#pivcol) { $c[$pivcol[$i]] = $M[$i][$nc]; }
  print "genus $h r $r k $k: ", ($consistent ? "IN SPAN, c_g = " . join(', ', map { $_->bstr } @c) : "NOT in span (" . scalar(@ks) . " monomials, $nc basis functions)"), "\n";
}
