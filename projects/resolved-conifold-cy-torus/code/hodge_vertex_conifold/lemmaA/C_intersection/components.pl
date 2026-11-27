#!/usr/bin/perl
# Extract the tau-independent Hodge components Xi^(k)_h(j_1..j_r), k = 0..h, from the exact table
# q_g<h>_r<r>.txt (Q = Sum_k s^k/k! <Omega prod lam_{j_i} Tt^k>, c-normalisation, s = -1 gives Xi_h).
# Method: at h+1 rational values of x put s = -1, rescale to c = 1 (factor c^(3h-3+2r)), and solve
#   Xi_h(j; tau) = Sum_{k=0}^h (-1/sigma)^k Xi^(k)(j),  sigma = -x/c^2,  c = 1 + x
# (brief (C)); one extra value of x is used as an independent consistency check.
# Output: the component polynomials (in j1..jr) written to comp_g<h>_r<r>.txt (one line per k).
# Usage: perl components.pl h r [datadir]
use strict; use warnings;
use FindBin; use lib $FindBin::Bin;
use polylib; use Math::BigRat;
my ($h, $r, $dir) = @ARGV; $dir //= "$FindBin::Bin/../data";
open(my $fh, '<', "$dir/q_g${h}_r${r}.txt") or die "no q_g${h}_r${r}.txt";
my $line = join('', <$fh>); close $fh;
my @f = split /\|/, $line; my $expr = $f[2]; $expr =~ s/\s+//g;
my @xs = map { Math::BigRat->new($_) } ('2', '3', '5', '7', '11', '13', '17', '19', '23', '29');
my $np = $h + 2; die "need more x values" if $np > @xs;
my (@t, @P);
for my $i (0 .. $np - 1) {
  my $x = $xs[$i]; my $c = $x + 1; my $sigma = -$x / ($c * $c);
  my $p = parse_expr($expr, { x => $x, s => -1 });
  $p = pscale($p, $c ** (3 * $h - 3 + 2 * $r));
  push @t, Math::BigRat->new(-1) / $sigma; push @P, $p;
}
# solve with the first h+1 points: Sum_k t_i^k X_k = P_i
my $n = $h + 1; my @A; my @B;
for my $i (0 .. $n - 1) { $A[$i] = [ map { $t[$i] ** $_ } 0 .. $h ]; $B[$i] = $P[$i]; }
for my $col (0 .. $n - 1) {
  my $piv = $col; $piv++ while $A[$piv][$col]->is_zero;
  @A[$col, $piv] = @A[$piv, $col]; @B[$col, $piv] = @B[$piv, $col];
  my $inv = Math::BigRat->new(1) / $A[$col][$col];
  $A[$col] = [ map { $_ * $inv } @{$A[$col]} ]; $B[$col] = pscale($B[$col], $inv);
  for my $row (0 .. $n - 1) { next if $row == $col; my $fct = $A[$row][$col]; next if $fct->is_zero;
    $A[$row] = [ map { $A[$row][$_] - $fct * $A[$col][$_] } 0 .. $n - 1 ]; $B[$row] = psub($B[$row], pscale($B[$col], $fct)); }
}
my @X = @B;
# consistency check at the extra point
my $chk = {}; for my $k (0 .. $h) { $chk = padd($chk, pscale($X[$k], $t[$n] ** $k)); }
my $diff = psub($chk, $P[$n]);
print STDERR "genus $h r $r: consistency check at extra x: ", (pis_zero($diff) ? "OK" : "FAIL"), "\n";
open(my $out, '>', "$FindBin::Bin/comp_g${h}_r${r}.txt") or die;
for my $k (0 .. $h) { print $out "$h | $r | $k | ", pstr($X[$k]), "\n"; }
close $out;
print STDERR "wrote comp_g${h}_r${r}.txt\n";
