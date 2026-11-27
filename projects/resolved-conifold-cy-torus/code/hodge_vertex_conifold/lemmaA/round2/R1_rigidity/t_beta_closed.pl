#!/usr/bin/perl
# t_beta_closed.pl -- compare the exact beta_{k,L}(tau) (polynomials in tau, from interp_beta.pl) with the closed form T1:
#   beta_{k,L} = kappa_k (-1)^d [ (1+tau) binom(N tau, d) - tau binom(N(1+tau)+d-1, d) ],  N = L+1, d = k-L,
#   kappa_1 = tau + 1/2, kappa_k = B_k/(k k!) (k >= 2).
# Polynomial identity in tau is checked at deg+2 distinct rational points (deg = max degree on either side).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my $file = shift // "$FindBin::Bin/out/beta_poly_K20.txt";
my @Bn = (R(1)); for my $m (1..40) { my $s = R(0); for my $k (0..$m-1) { $s = $s + binom($m+1, $k) * $Bn[$k]; } $Bn[$m] = -$s/($m+1); }
sub gbinom { my ($x, $d) = @_; my $v = R(1); for my $i (0..$d-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
open my $fh, '<', $file or die; my ($n, $bad, $kmax) = (0, 0, 0);
while (<$fh>) { next unless /^beta (\d+) (\d+) : (.*)$/; my ($k, $L) = ($1, $2); my @c = map { R($_) } split ' ', $3; $kmax = $k if $k > $kmax;
  my $deg = @c + ($k - $L) + 2;
  for my $i (0..$deg) { my $tau = R(2*$i + 3) / R(7*$i + 5) - R(1)/R(3);
    my $lhs = R(0); my $pw = R(1); for my $cc (@c) { $lhs = $lhs + $cc*$pw; $pw = $pw*$tau; }
    my $kap = ($k == 1) ? $tau + R('1/2') : $Bn[$k] / ($k * fact($k)); my $N = $L + 1; my $d = $k - $L;
    my $rhs = $kap * R((-1)**$d) * ((1+$tau)*gbinom($N*$tau, $d) - $tau*gbinom($N*(1+$tau)+$d-1, $d));
    $n++; $bad++ if $lhs != $rhs; } }
print "t_beta_closed: k <= $kmax, $n evaluations, $bad mismatches\n";
