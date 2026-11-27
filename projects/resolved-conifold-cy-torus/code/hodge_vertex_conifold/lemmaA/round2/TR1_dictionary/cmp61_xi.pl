#!/usr/bin/perl
# cmp61_xi.pl TAU SAVEFILE : compare the one-leaf TR values (xi-coefficients saved by run_xi.pl with SAVE=...),
#   Xi_h(j) = sigma^(h-1) Sum_k C_{h,1}[k] j^k      (P2, coefficient form; all j, including those with Pi_j(tau) = 0),
# with the conjectural closed formula (6.1) of proof_M_residue.md s.6:
#   sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j - l) [lambda^(2h)] S(e1 lambda)^(1 + j e2) S(e2 lambda)^(1 + j e1),
#   e1 = -tau, e2 = 1 + tau, S(y) = (y/2)/sinh(y/2),
# at j = 0..3h+1 (more than deg_j + 1 = 3h - 1 values, so agreement is an identity in j).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
my ($tau, $file) = (R(shift), shift);
my @B = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"), R("43867/798"), R("-174611/330"));
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 1..$n; $r }
my $sig = $tau * (1 + $tau);
sub pred { my ($h, $j) = @_; my $K = $h;
  my @lS = (R(0)); push @lS, -$B[$_] / (2 * $_ * fact(2 * $_)) for 1..$K;
  my $e1 = -$tau; my $e2 = 1 + $tau;
  my @c = (R(0)); push @c, (1 + $j * $e2) * $lS[$_] * ($e1 ** (2 * $_)) + (1 + $j * $e1) * $lS[$_] * ($e2 ** (2 * $_)) for 1..$K;
  my @r = (R(1)); push @r, R(0) for 1..$K;
  for my $n (1..$K) { my $s = R(0); $s += $_ * $c[$_] * $r[$n - $_] for 1..$n; $r[$n] = $s / $n; }
  my $p = R(1); $p *= ($j - $_) for 1..2*$h-2;
  return $p * $r[$h] / ($sig ** $h); }
open my $f, '<', $file or die "no $file";
while (my $line = <$f>) { next if $line =~ /^#/; chomp $line; my ($g, $n, $body) = split /\s*\|\s*/, $line, 3; next unless $n == 1 && $g >= 2;
  my %C; for my $t (split /\s*;\s*/, $body) { my ($k, $v) = split /:/, $t; $C{$k} = R($v); }
  my ($cnt, $bad, $zeros) = (0, 0, '');
  for my $jj (0..3*$g+1) { my $j = R($jj); my $tr = R(0); for my $k (keys %C) { $tr += $C{$k} * $j ** $k; } $tr *= $sig ** ($g - 1);
    my $p = pred($g, $j); $cnt++; $bad++ if $tr != $p; $zeros .= " $jj" if $tr->is_zero; }
  print "tau=$tau h=$g: compared Xi_h(j), j = 0..", 3*$g+1, " ($cnt values), mismatches with (6.1): $bad; zeros of Xi_h(j) at j =$zeros\n"; }
