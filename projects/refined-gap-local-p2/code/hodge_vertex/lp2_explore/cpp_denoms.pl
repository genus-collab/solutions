#!/usr/bin/perl
# cpp_denoms.pl -- list the distinct denominator factors of the identity-row coefficients cPP[0,0][k] stored in a cPP_z<K>.m file
# (Mathematica InputForm list). Purpose (CHECK 3): confirm that the only y-dependent denominator factor is
# (eps1 + eps2 + 3 LLL[0][y]) = E * 3t/(1+3t) (LLL -> LLt0 = -E/(3(1+3t))), i.e. that the identity row has no pole at any
# t != 0, in particular none at t = -1/3 (colour 0) or at t = -1/(3 w), -1/(3 w^2) (colours 1, 2 on the real branch).
# Usage: perl cpp_denoms.pl <cPP_z9.m>
use strict; use warnings;
local $/; my $s = <>;
$s =~ s/LLL\x5B0\x5D\x5By\x5D/L/g; $s =~ s/\x5C\x5BEpsilon\x5D\x5B(\d)\x5D/e$1/g; $s =~ s/\s+//g;
$s =~ s/^\{//; $s =~ s/\}$//;
# split at depth-0 commas
my @ent; my $d = 0; my $cur = '';
for my $c (split //, $s) { $d++ if $c eq '(' || $c eq '['; $d-- if $c eq ')' || $c eq ']';
  if ($c eq ',' && $d == 0) { push @ent, $cur; $cur = ''; } else { $cur .= $c; } }
push @ent, $cur;
print "entries: ", scalar(@ent), "\n";
my $yother = ($s =~ /y/) ? "yes" : "no"; print "y outside LLL[0][y]: $yother\n";
for my $k (0..$#ent) {
  my $e = $ent[$k]; my %fac;
  # every parenthesised group immediately after a '/'
  while ($e =~ m{/\(}g) { my $p = pos($e) - 1; my $dd = 0; my $j = $p;
    for (; $j < length $e; $j++) { my $c = substr($e,$j,1); $dd++ if $c eq '('; $dd-- if $c eq ')'; last if $dd == 0; }
    my $g = substr($e, $p+1, $j-$p-1);
    # depth-1 factors of the group separated by '*'
    my @f; my $dep = 0; my $cf = '';
    for my $c (split //, $g) { $dep++ if $c eq '(' || $c eq '['; $dep-- if $c eq ')' || $c eq ']';
      if ($c eq '*' && $dep == 0) { push @f, $cf; $cf = ''; } else { $cf .= $c; } }
    push @f, $cf;
    for (@f) { my $b = $_; $b =~ s/\^\d+$//; $fac{$b}++ if $b =~ /L/; $fac{"(L-free) $b"}++ if $b !~ /L/; } }
  my @Lf = grep { !/^\(L-free\)/ } sort keys %fac;
  my @nf = map { s/^\(L-free\) //r } grep { /^\(L-free\)/ } sort keys %fac;
  printf "k=%d: L-dependent denominator factors: %s | other: %s\n", $k, join(' ; ', @Lf), join(' ', @nf);
}
# degree check: max power of L in the numerator vs the power m of (e1+e2+3*L) in the denominator.
# With L = -E/(3(1+3t)) and E+3L = 3Et/(1+3t), deg_L N <= m means cPP[k] has no pole at t = -1/3 either (colour 0 finite).
for my $k (0..$#ent) {
  my $e = $ent[$k]; my $m = 0;
  my $den = $e; $m = $1 if $den =~ /\(e1\+e2\+3\*L\)\^(\d+)\)\s*$/ || $den =~ /\(e1\+e2\+3\*L\)\^(\d+)\)$/;
  $m = 1 if !$m && $e =~ /\(e1\+e2\+3\*L\)\)$/;
  (my $num = $e) =~ s/\(e1\+e2\+3\*L\)(\^\d+)?//g;
  my $dl = 0; while ($num =~ /L(\^(\d+))?/g) { my $p = defined $2 ? $2 : 1; $dl = $p if $p > $dl; }
  printf "k=%d: power of (E+3L) in denominator m = %d, max L-power elsewhere = %d -> %s\n", $k, $m, $dl, ($dl <= $m ? "finite at t=-1/3" : "POLE at t=-1/3");
}
