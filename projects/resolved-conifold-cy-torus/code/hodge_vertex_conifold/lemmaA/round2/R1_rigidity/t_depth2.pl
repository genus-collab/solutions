#!/usr/bin/perl
# t_depth2.pl -- closed form of the depth-2 closed-form coefficient beta_{L+2,L} (coefficient of y_L in relation_{L+2}):
#   beta_{L+2,L} = -(sigma/2)(L+1) [x^{L+2}](log(W/x) + x/W) - (sigma^2/2) [x^L] W'(x) (Int_0^x W)/W(x)^2,   W = e^{tau x}(e^x - 1),
# compared with the KP + leaf-curve computation (out/beta_tau2_7_K10.txt, tau = 2/7).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my $tau = R('2/7'); my $sig = $tau*(1+$tau); my $N = 14;
my @etx = map { $tau**$_ / fact($_) } 0..$N+1;            # e^{tau x}
my @em1 = (R(0), map { R(1)/fact($_) } 1..$N+1);           # e^x - 1
my $W = series_mul(\@etx, \@em1, $N+1);                     # W, W[0] = 0
my @Wx = map { $W->[$_+1] } 0..$N;                         # W/x
my $lWx = series_log(\@Wx, $N);
my $xW = series_inv(\@Wx, $N);                              # x/W
my @h = map { $lWx->[$_] + $xW->[$_] } 0..$N;
my @Wp = map { ($_+1) * $W->[$_+1] } 0..$N;                # W'
my @IW = (R(0), map { $W->[$_-1] / $_ } 1..$N+1);         # Int_0^x W  (starts at x^2)
# W'(x) Int W / W^2 = W' * (Int W/x^2) / (W/x)^2
my @IWx2 = map { $IW[$_+2] // R(0) } 0..$N;
my $q = series_mul(series_mul(\@Wp, \@IWx2, $N), series_inv(series_mul(\@Wx, \@Wx, $N), $N), $N);
my %tab; open my $fh, '<', "$FindBin::Bin/out/beta_tau2_7_K10.txt" or die; while (<$fh>) { if (/^REL (\d+): (.*)$/) { $tab{$1} = [ map { R($_) } split ' ', $2 ]; } } close $fh;
for my $L (1..8) { my $g = -($sig/2)*($L+1)*$h[$L+2] - ($sig*$sig/2)*$q->[$L];
  my $want = $tab{$L+2}[$L-1]; printf "L=%d  closed form %s  table %s  %s\n", $L, "$g", "$want", ($g == $want ? "OK" : "MISMATCH"); }
