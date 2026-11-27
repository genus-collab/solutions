#!/usr/bin/perl
# t_zexp.pl -- Perl twin of check (C4) of R1_checks.wls: T1 from an implementation of the genus-G linearised KP system along the
# leaf curve in the z-expansion (a = 1/z, b = 1/w), WITHOUT the y-coordinate Hankel form of ode_leaf.pl:
#   G(a,b) = a b QV(a,b) / (1 + a b sum_n c_n h_{n-1}(a,b)),  QV = sum_n V_n h_{n-1}(a,b),  h_{n-1} = sum_{i+j=n-1} a^i b^j,
#   c_n = sigma Pi_n W^(n+1)/(n+1)  (genus-0 curve at P(t): Y_P(z) = z - sigma sum_n Pi_n W^(n+1) z^-n/(n+1)),
#   V_n = n U_n + sum_j j P_j G_jn,  dU_n/dt = sum_i P_i' G_in,  dPhi/dt = sum_i P_i' U_i,  U_n(0) = y_{n-1}.
# Iterates the whole system K+2 times (each pass fixes one more t-order).  Exact (Math::BigRat).  Usage: perl t_zexp.pl K tau seed
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($K, $ts, $seed) = @ARGV; my $tau = R($ts); my $sig = $tau*(1+$tau); my $NM = $K + 2; srand($seed // 1);
sub gbin { my ($x, $m) = @_; return R(0) if $m < 0; my $v = R(1); for my $i (0..$m-1) { $v = $v * ($x - $i) / ($i + 1); } return $v; }
sub smul { series_mul($_[0], $_[1], $K) }
my $Tt = T_of_t_series($K); my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $K), $Tt, $K);
my (%Pj, %dPj);
for my $j (2..$NM) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $K), $Tt, $K); my $p = [R(1), map { R(0) } 1..$K]; $p = smul($p, $h1) for 1..$j; my $rho = smul($hj, series_inv($p, $K));
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$K ]; $dPj{$j} = [ map { ($_+1) * ($Pj{$j}[$_+1] // R(0)) } 0..$K ]; }
my @al = (R(1), map { R(0) } 1..$K);
for my $it (0..$K+1) { my @nw = (R(1), map { R(0) } 1..$K); my $wp = [@al]; for my $j (2..$NM) { $wp = smul($wp, \@al); my $q = smul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp); $nw[$_] = $nw[$_] + $q->[$_] for 0..$K; } @al = @nw; }
# trivariate truncated polynomials: hash "m,i,j" -> coeff (t^m a^i b^j), m <= K, i, j <= NM
sub tmul { my ($A, $B) = @_; my %C; for my $ka (keys %$A) { my ($m1, $i1, $j1) = split /,/, $ka; my $ca = $A->{$ka};
    for my $kb (keys %$B) { my ($m2, $i2, $j2) = split /,/, $kb; my ($m, $i, $j) = ($m1+$m2, $i1+$i2, $j1+$j2); next if $m > $K || $i > $NM || $j > $NM;
      my $k = "$m,$i,$j"; $C{$k} = ($C{$k} // R(0)) + $ca * $B->{$kb}; } }
  for (keys %C) { delete $C{$_} if $C{$_}->is_zero; } return \%C; }
sub tadd { my ($A, $B, $s) = @_; my %C = %$A; for my $k (keys %$B) { $C{$k} = ($C{$k} // R(0)) + ($s // R(1)) * $B->{$k}; } for (keys %C) { delete $C{$_} if $C{$_}->is_zero; } return \%C; }
# Emat = a b sum_n c_n(t) h_{n-1}(a,b)
my %E; my $wp = [@al];
for my $n (1..$NM) { $wp = smul($wp, \@al); my $cn = [ map { $sig * Pi_at($n, $tau) / ($n+1) * $_ } @$wp ];   # W^(n+1)
  for my $i (0..$n-1) { my $j = $n - 1 - $i; next if $i+1 > $NM || $j+1 > $NM; for my $m (0..$K) { next if $cn->[$m]->is_zero; my $k = "$m," . ($i+1) . "," . ($j+1); $E{$k} = ($E{$k} // R(0)) + $cn->[$m]; } } }
my $S = { "0,0,0" => R(1) }; for (1..$NM+1) { $S = tadd({ "0,0,0" => R(1) }, tmul(\%E, $S), R(-1)); }
my @y = (R(0), map { R(int(rand(19)) - 9) } 1..2*$NM+2);
my @U = map { [ map { R(0) } 0..$K ] } 0..$NM; for my $n (2..$NM) { $U[$n][0] = $y[$n-1]; }
my %G;   # "i,n" -> t-series
for my $pass (1..$K+2) {
  my @V; for my $n (1..$NM) { my $v = [ map { R($n) * $_ } @{$U[$n]} ]; for my $j (2..$NM) { my $g = $G{"$j,$n"}; next unless $g; my $q = smul($Pj{$j}, $g); $v->[$_] = $v->[$_] + R($j)*$q->[$_] for 0..$K; } $V[$n] = $v; }
  my %QV; for my $n (1..$NM) { for my $i (0..$n-1) { my $j = $n - 1 - $i; next if $i+1 > $NM || $j+1 > $NM; for my $m (0..$K) { next if $V[$n][$m]->is_zero; my $k = "$m," . ($i+1) . "," . ($j+1); $QV{$k} = ($QV{$k} // R(0)) + $V[$n][$m]; } } }
  my $Gab = tmul(\%QV, $S);
  %G = (); for my $k (keys %$Gab) { my ($m, $i, $j) = split /,/, $k; $G{"$i,$j"} //= [ map { R(0) } 0..$K ]; $G{"$i,$j"}[$m] = $Gab->{$k}; }
  for my $n (1..$NM) { my @d = map { R(0) } 0..$K; for my $i (2..$NM) { my $g = $G{"$i,$n"}; next unless $g; my $q = smul($dPj{$i}, $g); $d[$_] = $d[$_] + $q->[$_] for 0..$K; }
    my @u = map { R(0) } 0..$K; $u[0] = $n >= 2 ? $y[$n-1] : R(0); for my $m (0..$K-1) { $u[$m+1] = $d[$m] / ($m+1); } $U[$n] = \@u; } }
my @Phi = (R(0)); { my @d = map { R(0) } 0..$K; for my $i (2..$NM) { my $q = smul($dPj{$i}, $U[$i]); $d[$_] = $d[$_] + $q->[$_] for 0..$K; } for my $m (0..$K-1) { $Phi[$m+1] = $d[$m]/($m+1); } }
my @Bn = (R(1)); for my $mm (1..$K+2) { my $s = R(0); for my $k (0..$mm-1) { $s = $s + binom($mm+1, $k) * $Bn[$k]; } $Bn[$mm] = -$s/($mm+1); }
my $bad = 0;
for my $k (1..$K) { my $kap = $k == 1 ? $tau + R('1/2') : $Bn[$k]/($k*fact($k)); my $pred = R(0);
  for my $L (1..$k) { my $N = $L + 1; my $d = $k - $L; $pred = $pred + $y[$L] * R((-1)**$d) * ((1+$tau)*gbin($N*$tau, $d) - $tau*gbin($N*(1+$tau) + $d - 1, $d)); }
  $pred = $pred * $kap; $bad++ if $pred != $Phi[$k]; print "t^$k: z-expansion KP  $Phi[$k]   T1  $pred\n" if $ENV{VERBOSE} || $pred != $Phi[$k]; }
print "t_zexp K=$K tau=$ts: ", ($bad ? "$bad mismatches" : "T1 confirmed to order $K"), "\n";
