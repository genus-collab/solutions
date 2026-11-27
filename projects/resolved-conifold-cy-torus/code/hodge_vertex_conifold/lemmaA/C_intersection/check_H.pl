#!/usr/bin/perl
# Checks of Theorem H (Hurwitz component) and of its strong form H' on the component polynomials
# comp_g<h>_r<r>.txt produced by components.pl.
# The k-th component is X_k(j) = Sum_a <c(T(-log)) e_k(alpha+alpha^2) prod psi_i^a_i>_{h,r} prod j_i^a_i  (c = 1).
#  (1) lattice test: X_k(j) at all j_i >= 1 with D = Sum (j_i - 1) <= 2h-2 (range D <= 2h-3, plus the boundary D = 2h-2);
#  (2) strong form H': for every j in the range and every b with 0 <= b_i <= j_i - 1,
#        Int c(T(-log)) e_k prod_i psi_i^b_i / prod_{m=1}^{j_i} (1 - m psi_i)
#      = Sum_n prod_i h_{n_i}(1..j_i) I_k(b + n)   must vanish for k = 0 (Theorem H'); reported for every k.
# Usage: perl check_H.pl h r
use strict; use warnings;
use FindBin; use lib $FindBin::Bin;
use polylib; use Math::BigRat;
my ($h, $r) = @ARGV;
open(my $fh, '<', "$FindBin::Bin/comp_g${h}_r${r}.txt") or die; my @X;
while (<$fh>) { chomp; my @f = split /\s*\|\s*/; $X[$f[2]] = parse_expr($f[3], {}); } close $fh;
my $N = 3 * $h - 3 + $r; my @jv = map { "j$_" } 1 .. $r;
sub hcomp { my ($n, $j) = @_; my @H = (Math::BigRat->new(1), (Math::BigRat->new(0)) x $n); # h_m(1..0) = delta_m0
  for my $l (1 .. $j) { for my $m (1 .. $n) { $H[$m] = $H[$m] + $l * $H[$m - 1]; } } return $H[$n]; }
# enumerate j vectors (nondecreasing) with j_i >= 1 and D <= 2h-2
my @pts; my @cur;
sub rec { my ($i, $min) = @_; if ($i == $r) { my $D = 0; $D += $_ - 1 for @cur; push @pts, [@cur] if $D <= 2 * $h - 2; return; }
  for (my $j = $min; $j <= 2 * $h; $j++) { push @cur, $j; my $D = 0; $D += $_ - 1 for @cur; if ($D <= 2 * $h - 2) { rec($i + 1, $j); } pop @cur; } }
rec(0, 1);
print "genus $h, r = $r\n(1) lattice values X_k(j), j_i >= 1, D <= 2h-2\n";
for my $p (@pts) { my $D = 0; $D += $_ - 1 for @$p; my %v; @v{@jv} = @$p;
  my @vals = map { my $e = peval_at($X[$_], \%v); pcoef($e, '')->bstr } 0 .. $h;
  printf "  j = (%s)  D = %d %s  X_0..X_%d = %s\n", join(',', @$p), $D, ($D <= 2 * $h - 3 ? 'in range ' : 'BOUNDARY '), $h, join(' ; ', @vals); }
print "(2) strong form: Int c(T) e_k prod psi^b / prod_m (1 - m psi), all 0 <= b_i <= j_i - 1, j in range\n";
for my $k (0 .. $h) { my ($nz, $tot) = (0, 0); my @ex;
  for my $p (@pts) { my $D = 0; $D += $_ - 1 for @$p; next unless $D <= 2 * $h - 3;
    my @bs = ([]); for my $i (0 .. $r - 1) { my @nb; for my $b (@bs) { push @nb, [@$b, $_] for 0 .. $p->[$i] - 1; } @bs = @nb; }
    for my $b (@bs) { $tot++; my $S = Math::BigRat->new(0);
      for my $key (keys %{$X[$k]}) { my $e = mono_parse($key); my $ok = 1; my $w = Math::BigRat->new(1);
        for my $i (0 .. $r - 1) { my $a = $e->{$jv[$i]} // 0; my $nn = $a - $b->[$i]; if ($nn < 0) { $ok = 0; last; } $w *= hcomp($nn, $p->[$i]); }
        $S += $w * $X[$k]{$key} if $ok; }
      unless ($S->is_zero) { $nz++; push @ex, "j=(" . join(',', @$p) . ") b=(" . join(',', @$b) . "): " . $S->bstr if @ex < 6; } } }
  print "  k = $k: $nz non-zero out of $tot", (@ex ? "; e.g. " . join(' | ', @ex) : ''), "\n"; }
