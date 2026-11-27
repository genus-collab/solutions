#!/usr/bin/perl
# t_crit.pl -- consistency check: the constant C(t) of the curve identity u_t(phi_t(lambda)) = u(lambda)/C equals
#   C_a = alpha e^{tau^2 t}/H_1(T)   (Lagrange inversion + Pfaff)   and   C_* = u(1+tau)/u_t(mu_c)   (critical values),
# mu_c the critical point of u_t:  e^{t mu_c} = Lc,  (1+tau) Lc^2 - (1+2 tau) E Lc + tau E = 0,  E = e^t,
#   Lc = E (1 + 2 tau + S)/(2(1+tau)),  S = sqrt(1 + 4 sigma T),  T = 1 - e^-t.
# log u_t(mu_c) = log t - sigma t + tau log(Lc - E) - (1+tau) log(Lc - 1) - tau log Lc.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my ($taus, $M) = @ARGV; my $tau = R($taus); my $sig = $tau*(1+$tau); my $N = $M + 3;
my @E = map { R(1)/fact($_) } 0..$N; my $Tt = T_of_t_series($N);
my @a = map { 4*$sig*$Tt->[$_] } 0..$N; $a[0] = R(1); my $S = series_pow(\@a, R('1/2'), $N);
my @num = map { $S->[$_] } 0..$N; $num[0] = $num[0] + 1 + 2*$tau; my $Lc = series_mul(\@E, [ map { $_ / (2*(1+$tau)) } @num ], $N);
my @LmE = map { $Lc->[$_] - $E[$_] } 0..$N; my @Lm1 = @$Lc; $Lm1[0] = $Lm1[0] - 1;
die "Lc-E not O(t)" unless $LmE[0]->is_zero && $Lm1[0]->is_zero;
my @q1 = map { $LmE[$_+1] / $tau } 0..$N-1; my @q2 = map { $Lm1[$_+1] / (1+$tau) } 0..$N-1;   # (Lc-E)/(tau t), (Lc-1)/((1+tau) t)
my $l1 = series_log(\@q1, $M); my $l2 = series_log(\@q2, $M); my $l3 = series_log($Lc, $M);
# log C_* = log u(1+tau) - log u_t(mu_c) = -[ -sigma t + tau l1 - (1+tau) l2 - tau l3 ]  (the t, tau^tau, (1+tau)^(-1-tau) factors cancel)
my @lCs = map { -($tau*$l1->[$_] - (1+$tau)*$l2->[$_] - $tau*$l3->[$_]) } 0..$M; $lCs[1] = $lCs[1] + $sig;
# C_a from the leaf curve
my $h1 = series_compose(hyp2f1_series(-$tau, -$tau, R(1), $N), $Tt, $N);
my %Pj; for my $j (2..$N+1) { my $hj = series_compose(hyp2f1_series(-$j*$tau, -$j*$tau, R($j), $N), $Tt, $N); my $p = [R(1), map { R(0) } 1..$N]; $p = series_mul($p, $h1, $N) for 1..$j; my $rho = series_mul($hj, series_inv($p, $N), $N);
  my $cc = R((-1)**$j) * Pi_at($j, $tau) / $j; $Pj{$j} = [ map { $_ >= $j-1 ? $cc * $rho->[$_-$j+1] : R(0) } 0..$N ]; }
my @al = (R(1), map { R(0) } 1..$N);
for my $it (0..$N+1) { my @nw = (R(1), map { R(0) } 1..$N); my $wp = [@al]; for my $j (2..$N+1) { $wp = series_mul($wp, \@al, $N); my $q = series_mul([ map { Pi_at($j,$tau) * $_ } @{$Pj{$j}} ], $wp, $N); $nw[$_] = $nw[$_] + $q->[$_] for 0..$N; } @al = @nw; }
my $lCa = series_log(series_mul(\@al, series_inv($h1, $N), $N), $M); $lCa->[1] = $lCa->[1] + $tau*$tau;
my $ok = 1; for my $m (0..$M) { $ok = 0 if $lCs[$m] != $lCa->[$m]; }
print "tau=$taus: log C_* = ", join(' ', map { "$_" } @lCs[0..$M]), "\n         log C_a = ", join(' ', map { "$_" } @$lCa[0..$M]), "  ", ($ok ? "EQUAL" : "DIFFERENT"), "\n";
