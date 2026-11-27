#!/usr/bin/perl
# pk_fit.pl -- structure of the normalised closed-form relations in the F-basis:
#   (k+1) R_k = [y^{k+1}] F(y) (1 + tau A_k(y) + tau^2 B_k(y)),   check A_k - B_k = 1 + y^k - (1-y)^k,
#   B_k(y) = -(k+1) y^2 P_k(y), P_k palindromic of degree k-2, P_k = y^m p_k(v), v = (1-y)^2/y = y + 1/y - 2, m = (k-2)/2.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my $file = shift // "$FindBin::Bin/out/beta_poly_K14.txt"; my $Kmax = shift // 14;
my @lines = `perl "$FindBin::Bin/R_bases.pl" "$file" $Kmax`;
for my $ln (@lines) { next unless $ln =~ /^R_(\d+) \(F-basis\): (.*)$/; my $k = $1; my $rest = $2; my (%A, %B);
  while ($rest =~ /F(\d+):(tau-poly\[[^\]]*\]|\[[^\]]*\])/g) { my ($n, $p) = ($1, $2); my $m = $k + 1 - $n; my @c;
    if ($p =~ /^tau-poly\[(.*)\]$/) { @c = map { R($_) } split /, /, $1; } else { $p =~ /^\[(.*)\]$/; @c = (R($1)); }
    $A{$m} = ($c[1] // R(0)) * ($k+1); $B{$m} = ($c[2] // R(0)) * ($k+1); }
  my $ok = 1; for my $m (1..$k-1) { my $want = R((-1)**($m+1)) * binom($k, $m); $ok = 0 if ($A{$m} // R(0)) - ($B{$m} // R(0)) != $want; }
  # P_k coefficients: B_k = -(k+1) y^2 P_k  -> P_k[j] = -B[j+2]/(k+1), j = 0..k-3 (j = k-2 from palindromy)
  my $mm = ($k - 2)/2; my @P = map { -($B{$_+2} // R(0))/($k+1) } 0..$k-3; push @P, $P[0];
  # p_k(v): P_k(y) = y^mm p_k(v) with v = (1-y)^2/y: write p_k(v) = sum_j e_j v^j, i.e. P_k = sum_j e_j (1-y)^{2j} y^{mm-j}; solve triangularly from low powers of y
  my @e; my @res = @P;   # residual polynomial in y
  for my $j (reverse 0..$mm) { # coefficient of y^{mm-j} lowest power from term j=mm is y^0: process j = mm first (lowest y power 0)
  }
  # solve: for j from mm down to 0, the lowest remaining power is y^{mm-j}
  for (my $j = $mm; $j >= 0; $j--) { my $pw = $mm - $j; my $c = $res[$pw]; $e[$j] = $c;
    # subtract c (1-y)^{2j} y^{pw}
    for my $i (0..2*$j) { $res[$pw+$i] = $res[$pw+$i] - $c * binom(2*$j, $i) * R((-1)**$i); } }
  my $resid = grep { !$_->is_zero } @res;
  printf "k=%2d: A-B binomial check %s; p_k(v) = %s   (residual %d)\n", $k, ($ok ? "OK" : "FAIL"), join(' + ', map { "($e[$_])v^$_" } 0..$mm), $resid;
}
