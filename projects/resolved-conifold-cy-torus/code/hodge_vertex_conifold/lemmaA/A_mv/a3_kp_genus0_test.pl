#!/usr/bin/perl
# a3_kp_genus0_test.pl -- validation of the KP-relation generator of a3_kp_nullity.pl on genus 0, where everything is explicit:
# the genus-0 dilaton-point values Xi_0(T) = Gamma(J)/Gamma(J-r+3) (unit parts included) must satisfy the dispersionless Hirota equation
#   [z^-al w^-be p_R] exp( sum_{a,b} z^-a w^-b om0_ab ) = om0_{1,al+be-1}(R),   om0_ab(R) = sigma Pi_a Pi_b Xi_0(R + {a,b}),
# for every R and (al, be), exactly (all factors genus 0).  Exact rational arithmetic at a given tau.
# Usage: perl a3_kp_genus0_test.pl tau maxw maxR
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($tau, $maxw, $maxR) = @ARGV; $tau = R($tau); my $sig = $tau*(1+$tau);
sub Pi { my $j = shift; return R(1) if $j <= 1; my $p = R(1); for my $b (1..$j-1) { $p *= ($j*$tau + $b); } return $p / fact($j-1); }
sub X0 { my @L = @_; my $J = 0; $J += $_ for @L; my $r = scalar @L; if ($r >= 3) { my $v = R(1); $v *= ($J - $_) for 1..($r-3); return $v; } elsif ($r == 2) { return R(1)/$J; } else { return R(1)/($J*($J+1)); } }
my @Rs = ([]); for my $a (2..$maxR) { push @Rs, [$a]; for my $b (2..$a) { push @Rs, [$a,$b]; } }
my ($n, $bad) = (0, 0);
for my $R (@Rs) { my @R = @$R; my $nR = scalar @R; my $full = (1 << $nR) - 1;
  my %om; for my $B (0..$full) { my @RB = map { $R[$_] } grep { ($B >> $_) & 1 } 0..$nR-1; for my $a (1..$maxw) { for my $b (1..$maxw) { $om{"$a,$b,$B"} = $sig * Pi($a) * Pi($b) * X0(@RB, $a, $b); } } }
  my %Em; my $E; $E = sub { my ($a, $b, $S) = @_; my $k = "$a,$b,$S"; return $Em{$k} if exists $Em{$k}; my $res = R(0);
    if ($a == 0 && $b == 0) { $res = $S == 0 ? R(1) : R(0); } elsif ($a <= 0 || $b <= 0) { $res = R(0); }
    else { for my $a1 (1..$a) { for my $b1 (1..$b) { my $B = $S; while (1) { my $e = $E->($a-$a1, $b-$b1, $S & ~$B); $res += $a1 * $om{"$a1,$b1,$B"} * $e unless $e->is_zero; last if $B == 0; $B = ($B - 1) & $S; } } } $res /= $a; }
    $Em{$k} = $res; return $res; };
  for my $al (1..$maxw) { for my $be ($al..$maxw) { next if $al + $be > $maxw + 1;
      my $lhs = $E->($al, $be, $full); my $c = $al + $be - 1; my $rhs = $sig * Pi($c) * X0(@R, 1, $c); $n++;
      if ($lhs != $rhs) { $bad++; print "MISMATCH R=(@R) al=$al be=$be: lhs=$lhs rhs=$rhs\n" if $bad <= 5; } } } }
print "genus-0 dilaton-point dispersionless Hirota: $n relations tested (tau=$tau), $bad mismatches\n";
