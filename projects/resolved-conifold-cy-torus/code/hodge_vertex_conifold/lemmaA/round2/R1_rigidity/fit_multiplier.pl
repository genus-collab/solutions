#!/usr/bin/perl
# fit_multiplier.pl -- test whether the state along the leaf curve is a Toeplitz multiplication of the initial datum:
#   Uhat_t(lambda) = [ M_t(lambda) Psi(lambda) ]_-   (negative powers of lambda),  M_t universal (independent of the solution).
# For several basic solutions Psi = u^N (N = L0+1), fit M_m(lambda) (polynomial in lambda) at each t-order m from the exact state
# (ode_state.pl) using the valid coefficients (level <= K+1), then compare the fits across N.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($K, $ts, @L0s) = @ARGV; @L0s = (2, 3, 4) unless @L0s; my $tau = R($ts);
my %fits;
for my $L0 (@L0s) {
  my @lines = `perl "$FindBin::Bin/ode_state.pl" $K $ts $L0`;
  my %st; for (@lines) { next unless /^t\^(\d+) Uhat\(y\): (.*)$/; my $m = $1; for my $e (split /\s+/, $2) { next unless $e =~ /^y(\d+):(.*)$/; $st{"$m,$1"} = R($2); } }
  # Psi = u^N with u = y(1-y)^tau, as series in y (= 1/lambda)
  my $N = $L0 + 1; my $NN = $K + 4; my $u = series_mul([R(0), R(1)], series_pow([R(1), R(-1)], $tau, $NN), $NN); my $psi = [R(1), map { R(0) } 1..$NN]; $psi = series_mul($psi, $u, $NN) for 1..$N;
  for my $m (1..$K) {
    # unknown M_m(lambda) = sum_{p=0}^{m} a_p lambda^p ; [lambda^p * Psi]_- coefficient at y^n is psi[n+p] (n >= 1)
    # valid data: y^n with (n-1)+m <= K+1  and n >= 1
    my @rows; for my $n (1..$K+2-$m) { my @r = map { $psi->[$n+$_] // R(0) } 0..$m; push @r, $st{"$m,$n"} // R(0); push @rows, \@r; }
    # least-squares not needed: solve exactly (Gaussian elimination) and check consistency
    my $nc = $m + 1; my @A = map { [ @$_ ] } @rows; my $rk = 0; my @piv;
    for my $c (0..$nc-1) { my $pr; for my $i ($rk..$#A) { unless ($A[$i][$c]->is_zero) { $pr = $i; last; } } next unless defined $pr; @A[$rk, $pr] = @A[$pr, $rk];
      my $inv = R(1)/$A[$rk][$c]; $A[$rk] = [ map { $_ * $inv } @{$A[$rk]} ]; for my $i (0..$#A) { next if $i == $rk || $A[$i][$c]->is_zero; my $f = $A[$i][$c]; $A[$i] = [ map { $A[$i][$_] - $f*$A[$rk][$_] } 0..$nc ]; } push @piv, $c; $rk++; }
    my $incons = 0; for my $i ($rk..$#A) { $incons++ unless $A[$i][$nc]->is_zero; }
    my @sol = map { R(0) } 0..$m; for my $i (0..$rk-1) { $sol[$piv[$i]] = $A[$i][$nc]; }
    $fits{"$L0,$m"} = [ \@sol, $incons, $rk ];
    printf "L0=%d t^%d: M_m = %s   (rank %d of %d, inconsistent rows %d)\n", $L0, $m, join(' + ', map { "($sol[$_])l^$_" } 0..$m), $rk, $nc, $incons; } }
