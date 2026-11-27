#!/usr/bin/perl
# t_Rformula.pl -- test of the closed form of the closed-form (polar-order) relations:
#   relation_k = kappa_k R_k,   sum_k (k+1) R_k s^{k+1} = (1 + tau - sigma s/(1-s)) F(s) + (sigma s/(1+s) - tau) F(s/(1+s)),
# F(y) = sum_c c y_{c-1} u(y)^c, u = y(1-y)^tau  (y_L = Pi_{L+1} Xi(L+1)).  Compared with the exact rows beta_{k,L}(tau) (out/beta_poly_K14.txt).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my %B; open my $fh, '<', "$FindBin::Bin/out/beta_poly_K14.txt" or die; while (<$fh>) { next unless /^beta (\d+) (\d+) : (.*)$/; $B{"$1,$2"} = [ map { R($_) } split ' ', $3 ]; } close $fh;
my $K = 14; my $N = $K + 2;
my @Bn = (R(1)); for my $m (1..40) { my $s = R(0); for my $k (0..$m-1) { $s = $s + binom($m+1, $k) * $Bn[$k]; } $Bn[$m] = -$s/($m+1); }
my ($n, $bad) = (0, 0);
for my $taus ('2/7', '-3/11', '5/2', '1/3') { my $tau = R($taus); my $sig = $tau*(1+$tau);
  srand(17); my @y = (R(0), map { R(int(rand(21)) - 10)/R(1 + int(rand(5))) } 1..$K);   # y_1..y_K random
  # beta rows evaluated at tau
  my %beta; for my $k (1..$K) { for my $L (1..$k) { my $p = $B{"$k,$L"} // [R(0)]; my $v = R(0); my $pw = R(1); for my $c (@$p) { $v = $v + $c*$pw; $pw = $pw*$tau; } $beta{"$k,$L"} = $v; } }
  # F(y) = sum_c c y_{c-1} u^c, u = y(1-y)^tau
  my $u = series_mul([R(0), R(1)], series_pow([R(1), R(-1)], $tau, $N), $N);
  my @F = map { R(0) } 0..$N; my $up = [R(1), map { R(0) } 1..$N];
  for my $c (1..$K+1) { $up = series_mul($up, $u, $N); next if $c < 2; my $w = R($c) * $y[$c-1]; for my $i (0..$N) { $F[$i] = $F[$i] + $w * $up->[$i]; } }
  # R(s) = (1 + tau - sigma s/(1-s)) F(s) + (sigma s/(1+s) - tau) F(s/(1+s))
  my @k1 = (1 + $tau, map { -$sig } 1..$N);                       # 1 + tau - sigma sum_{j>=1} s^j
  my $t1 = series_mul(\@k1, \@F, $N);
  my @ms = (R(0), map { R((-1)**($_+1)) } 1..$N);                  # s/(1+s)
  my $Fc = series_compose(\@F, \@ms, $N);
  my @k2 = (-$tau, map { $sig * R((-1)**($_+1)) } 1..$N);         # sigma s/(1+s) - tau
  my $t2 = series_mul(\@k2, $Fc, $N);
  for my $k (1..$K) { my $pred = ($t1->[$k+1] + $t2->[$k+1]) / ($k+1);
    my $rel = R(0); for my $L (1..$k) { $rel = $rel + $beta{"$k,$L"} * $y[$L]; }
    my $kap = $k == 1 ? $tau + R('1/2') : $Bn[$k]/($k*fact($k));
    $n++; my $ok = ($k % 2 == 0 || $k == 1) ? ($rel == $kap * $pred) : ($rel->is_zero);
    if ($k == 1) { $ok = ($rel == $kap * $y[1]); }
    unless ($ok) { $bad++; print "tau=$taus k=$k: relation $rel, kappa*pred ", $kap*$pred, "\n"; } }
}
print "t_Rformula: $n comparisons, $bad mismatches\n";
