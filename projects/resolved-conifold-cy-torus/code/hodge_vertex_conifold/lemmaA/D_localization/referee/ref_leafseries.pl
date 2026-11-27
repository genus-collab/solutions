# ref_leafseries.pl -- the closed-form leaf series of leaf_structure_proof.md Theorem 11.1 (same construction as ref_closedform.pl):
#   leaf_series($x, $N) -> (\%A, \%C): A{j} = t-series (orders 0..N) of the lam_j coefficient of rho q_t, C{j} = C_j
#   G_R(\%A, $c, @R)    -> t-series of the coefficient of the irreducible symbol Xi^c(R) in F_h(rho q_t) (lam_1 parts resummed by the
#                          string identity Xi^c(1^p, R) = (J_R)_p c^(-2p) Xi^c(R)); R = () gives 1 (Theorem C: Xi(1^p) = 0 for p >= 1)
use strict; use warnings; use FindBin; require "$FindBin::Bin/ref_lib.pl";
sub leaf_series { my ($x, $N) = @_; my $M = $N + 2; my $c = 1 + $x;
  my $tau = s_zero($M); $tau->[$_] = rq((-1)**($_+1)) / rfact($_) for 1..$M;
  my $mu = s_scale($tau, -4 * $x / $c**2);
  my $sq = s_pow_rat($mu, rq(1,2));
  my $sstar = s_scale(s_inv(s_add(s_const(1, $M), $sq)), 2 / $c);
  my $dphi = s_int(s_sub(s_const(1, $M), s_scale($sstar, $x)));
  my $zt = s_shift_down(s_mul($tau, s_exp(s_scale($dphi, 1 / $c))));
  my @taupow = (s_const(1, $M)); push @taupow, s_mul($taupow[-1], $tau) for 1..$M;
  my (%A, %C);
  for my $j (1..$N+1) { my $P = rq(1); $P *= ($j - $_ * $c)**2 for 1..$j-1;
    $C{$j} = $P / (rfact($j) * rfact($j-1)) * $c**(4 - 2*$j);
    my $F = s_zero($M); my $a = $j * $x / $c;
    for my $n (0..$M) { $F = s_add($F, s_scale($taupow[$n], poch($a, $n)**2 / (poch(rq($j), $n) * rfact($n)))); }
    my $core = s_scale(s_mul(s_powint($zt, $j), $F), rq((-1)**$j) * $C{$j});
    my $Aj = s_zero($N); for my $k ($j-1..$N) { $Aj->[$k] = $core->[$k - $j + 1]; } $A{$j} = $Aj; }
  return (\%A, \%C); }
sub G_R { my ($A, $c, @R) = @_; my $N = $#{$A->{1}}; my $J = 0; $J += $_ for @R;
  my $Y1 = [ @{$A->{1}} ]; $Y1->[0] += $c**2;
  my $G = s_pow_rat(s_scale($Y1, -1 / $c**2), rq(-$J));
  my %m; $m{$_}++ for @R; my $aut = rq(1);
  for my $j (keys %m) { die "G_R: part $j < 2" if $j < 2; $G = s_mul($G, s_powint($A->{$j}, $m{$j})); $aut *= rfact($m{$j}); }
  return s_scale($G, 1 / $aut); }
1;
