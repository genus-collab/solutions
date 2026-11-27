#!/usr/bin/perl
# Check of the strong form H' of Theorem H (k = 0 only):
#   Int c_m(-R pi_* S_D) prod psi_i^{b_i} = 0  for all 0 <= b_i <= j_i with Sum b_i <= J - 1  (m = N - Sum b_i >= N - J + 1),
# where c(-R pi_* S_D) = c(T(-log)) / prod_i prod_{k=1}^{j_i} (1 - k psi_i). Uses the k = 0 component polynomial
# X_0(j) = Sum_a <c(T(-log)) prod psi^a> prod j^a from comp_g<h>_r<r>.txt. Also reports the excluded case b = j (all i), which is the
# honest top Chern class of the cokernel and need not vanish.  Usage: perl check_Hprime.pl h r
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use polylib; use Math::BigRat;
my ($h, $r) = @ARGV;
open(my $fh, '<', "$FindBin::Bin/comp_g${h}_r${r}.txt") or die; my $X;
while (<$fh>) { chomp; my @f = split /\s*\|\s*/; $X = parse_expr($f[3], {}) if $f[2] == 0; } close $fh;
my @jv = map { "j$_" } 1 .. $r;
sub hcomp { my ($n, $j) = @_; my @H = (Math::BigRat->new(1), (Math::BigRat->new(0)) x $n);
  for my $l (1 .. $j) { for my $m (1 .. $n) { $H[$m] = $H[$m] + $l * $H[$m - 1]; } } return $H[$n]; }
my @pts; my @cur;
sub rec { my ($i, $min) = @_; if ($i == $r) { my $D = 0; $D += $_ - 1 for @cur; my $J = 0; $J += $_ for @cur; push @pts, [@cur] if ($D <= 2 * $h - 3 && $J >= 1); return; }
  for (my $j = $min; $j <= 2 * $h; $j++) { push @cur, $j; rec($i + 1, $j); pop @cur; } }
rec(0, 0);
my ($tot, $nz, $exc) = (0, 0, 0); my @ex;
for my $p (@pts) { my $J = 0; $J += $_ for @$p;
  my @bs = ([]); for my $i (0 .. $r - 1) { my @nb; for my $b (@bs) { push @nb, [@$b, $_] for 0 .. $p->[$i]; } @bs = @nb; }
  for my $b (@bs) { my $sb = 0; $sb += $_ for @$b; my $S = Math::BigRat->new(0);
    for my $key (keys %$X) { my $e = mono_parse($key); my $ok = 1; my $w = Math::BigRat->new(1);
      for my $i (0 .. $r - 1) { my $a = $e->{$jv[$i]} // 0; my $nn = $a - $b->[$i]; if ($nn < 0) { $ok = 0; last; } $w *= hcomp($nn, $p->[$i]); }
      $S += $w * $X->{$key} if $ok; }
    if ($sb <= $J - 1) { $tot++; unless ($S->is_zero) { $nz++; push @ex, "j=(@$p) b=(@$b): " . $S->bstr if @ex < 5; } }
    else { $exc++; } } }
print "genus $h r $r: H' tested on $tot (j,b) pairs (j_i >= 0, J >= 1, D <= 2h-3, b_i <= j_i, Sum b <= J-1): $nz non-zero",
  (@ex ? " e.g. " . join(' | ', @ex) : ''), "; excluded b = j cases: $exc\n";
