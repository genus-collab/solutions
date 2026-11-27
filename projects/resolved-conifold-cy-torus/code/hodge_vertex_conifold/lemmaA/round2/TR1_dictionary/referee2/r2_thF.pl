#!/usr/bin/perl
# r2_thF.pl TAU E G : one-point test of Theorem F (the TR family y_e = log(1+w) - e/w) at genus G, against the Route E referee's
# TR-free engine (RefHodge::xi: dilaton leaves as explicit marked points, no kappa classes).
#   TR side : TRxi.pm (author's module, copied unchanged) with y = 'mix', e = E; one-point polynomial P_TR(j) = Sum_k C_{G,1}[k] j^k.
#   Hodge   : q[m](j) = (1/m!) Int_{Mbar_{G,1+m}} Omega lam_j(psi_1) prod Tt(psi),  Tt = psi^2/(1-psi);
#             Theorem F predicts  P_TR(j) = sigma^G (1 + E sigma)^(1-2G) Sum_m (-s)^m q[m](j),   s = E sigma/(1 + E sigma).
#   Also reported from the same Hodge data: E = 0 (BKMP: P(j) = sigma^G q[0](j)) is the s = 0 value, and the dilaton point
#   s = 1 (P2) is Sum_m (-1)^m q[m](j).  Compared at j = 0..3G+1 (> deg_j + 1 = 3G - 1 values: an identity in j).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Time::HiRes qw(time);
use TRxi; use RefHodge;
sub R { Math::BigRat->new(@_) }
$| = 1; my ($tau, $E, $G) = @ARGV; my $t0 = time;
my $tauR = R($tau); my $eR = R($E); my $sig = $tauR * (1 + $tauR);
my @J = map { R($_) } 0..3*$G+1;
my %res;
for my $case ('mix', 'bkmp', 'base') {
  my $kmax = 3*$G - 2; TRxi::setup(tau => $tau, kmax => $kmax, y => $case, e => $E);
  TRxi::run_all($G + 1, $G); my $W = $TRxi::W{"$G,1"};
  $res{$case} = [ map { my $j = $_; my $s = R(0); for my $k (keys %$W) { $s += $W->{$k} * $j ** $k; } $s } @J ];
  printf "TR %s done (%.1f s)\n", $case, time - $t0; }
my $om = RefHodge::omega($G, [R(1), $tauR, -1 - $tauR], 3*$G - 3);
my @Q = map { RefHodge::xi($G, $om, [$_], undef) } (0..3*$G+1);
printf "Hodge side (RefHodge::xi) done (%.1f s)\n", time - $t0;
my $s = $eR * $sig / (1 + $eR * $sig);
my ($bm, $bb, $bp) = (0, 0, 0);
for my $i (0..$#J) { my $q = $Q[$i];
  my $vm = R(0); my $pw = R(1); for my $m (0..$#$q) { $vm += $q->[$m] * $pw; $pw *= -$s; }
  my $predm = $sig ** $G * (1 + $eR * $sig) ** (1 - 2*$G) * $vm;
  my $predb = $sig ** $G * $q->[0];
  my $vp = R(0); $pw = R(1); for my $m (0..$#$q) { $vp += $q->[$m] * $pw; $pw *= R(-1); }
  my $predp = $sig ** (1 - $G) * $vp;
  $bm++ if $res{mix}[$i] != $predm; $bb++ if $res{bkmp}[$i] != $predb; $bp++ if $res{base}[$i] != $predp; }
printf "tau = %s, sigma = %s, genus %d, one leaf, j = 0..%d:\n  Theorem F at e = %s (s = %s): %d mismatches\n  BKMP (e = 0): %d mismatches\n  P2 (base curve): %d mismatches\n  time %.1f s\n",
  $tauR, $sig, $G, 3*$G+1, $eR, $s, $bm, $bb, $bp, time - $t0;
