#!/usr/bin/perl
# timing of the set-up phases of TRtau.pm (diagnostic)
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use Math::BigRat; use Time::HiRes qw(time);
use TRtau;
sub R { Math::BigRat->new(@_) }
my ($tau, $L) = @ARGV; $tau //= 1; $L //= 32;
my $t0 = time; my $tauR = R($tau);
my @X; my $t1 = $tauR / ($tauR + 1);
for my $j (2..$L+3) { $X[$j] = R((-1) ** ($j - 1)) / $j * ($tauR ** $j - ($tauR + 1) * $t1 ** $j); }
my $c2 = -$X[2]; my @q = map { -$X[$_ + 2] / $c2 } 0..$L+1;
printf "X coefficients: %.2f s\n", time - $t0; $t0 = time;
my @u = (R(0)); for my $n (1..$L+1) { my $f = TRtau::ps_powr(\@q, R(-$n) / 2, $n - 1); $u[$n] = $f->[$n-1] / $n; }
printf "Lagrange inversion: %.2f s; u_10 = %s (length of string %d)\n", time - $t0, $u[10]->bstr, length($u[$L]->bstr); $t0 = time;
my @uz = @u[1..$#u]; my $r = TRtau::ps_inv(\@uz, $L);
printf "inverse: %.2f s\n", time - $t0; $t0 = time;
my $p = TRtau::ps_mul($r, $r, $L);
printf "one product: %.2f s\n", time - $t0;
