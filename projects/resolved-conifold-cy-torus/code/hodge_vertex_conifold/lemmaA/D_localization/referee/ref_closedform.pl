# ref_closedform.pl -- referee check of the leaf coefficients behind Route D section 5 (Lemma 5.2), from the CLOSED FORM of
# leaf_structure_proof.md Theorem 11.1, computed here from scratch as exact t-series at rational x:
#   tau = 1 - e^-t, u = 4 x tau/c^2, s* = 2/(c (1 + (1-u)^(1/2))), dphi' = 1 - x s* (dphi(0) = 0), zeta = tau e^(dphi/c),
#   A_j(t) = (-1)^j C_j t^-1 zeta^j 2F1(jx/c, jx/c; j; tau),  C_j = P_{j-1}/(j! (j-1)! c^(2j-4)),  P_k = prod_{i=1}^k (k+1-ic)^2,
#   rho = (tau/t) 2/((1-u)^(1/4) (1 + (1-u)^(1/2))).
# Checks: (1) sum rule; (2) Sum j A_j = -c^2 rho; (3) rho and all a_j^(k), k <= 5, against the Mathematica fits
# hvc_leaf_structure_data_kt5.m (independent data); (4) top / sub-top coefficients (Cor. 11.3) to order N; (5) C_j = (c^2/j) eps_j(sigma);
# (6) Lemma 5.2 in resummed form: for every irreducible R (parts >= 2) with D(R) <= N - 1 the coefficient series of Xi(R) in F_h(rho q_t),
#     G_R(t) = (1 - Y_1/c^2)^(-J_R) prod_{j>=2} Y_j^(m_j)/m_j!   (Y_1 = A_1 + c^2, Y_j = A_j; the factor resums the lam_1 parts by the
#     string identity Xi(1^p, R) = (J_R)_p c^(-2p) Xi(R), genus independent), satisfies [t^D] G_R = prod (-1)^j C_j/|Aut R| and [t^(D+1)] G_R = 0;
#     [t^(D+2)] G_R is printed (non-zero = the deficit-2 terms that Theorem D5 kills by induction).
use strict; use warnings; use FindBin; require "$FindBin::Bin/ref_lib.pl";
my $N = 8;           # t-order of all checks
my $M = $N + 2;      # working order (margin for the division by t)
my $kt5file = "$FindBin::Bin/../../../hvc_leaf_structure_data_kt5.m";
my $kt5raw = do { local $/; open my $fh, '<', $kt5file or die "$kt5file: $!"; <$fh> };
# keep only the first two top-level elements {rho coefficients, dilaton decompositions} (the rest has rules and Missing[...])
sub top_elements { my ($s) = @_; $s =~ s/\s+//g; die "not a list" unless substr($s, 0, 1) eq '{';
  my ($depth, $start, @el) = (0, 1);
  for my $i (1..length($s)-1) { my $ch = substr($s, $i, 1);
    if ($ch eq '{' || $ch eq '(') { $depth++ } elsif ($ch eq '}' || $ch eq ')') { if ($depth == 0) { push @el, substr($s, $start, $i - $start); last; } $depth-- }
    elsif ($ch eq ',' && $depth == 0) { push @el, substr($s, $start, $i - $start); $start = $i + 1; } }
  return @el; }
my @kt5el = top_elements($kt5raw);
my $kt5src = "{" . $kt5el[0] . "," . $kt5el[1] . "}";
my $bad = 0;
sub ok { my ($name, $a, $b) = @_; my $g = ($a == $b); $bad++ unless $g; print "   $name: ", ($g ? "ok" : "MISMATCH ($a vs $b)"), "\n" unless $g && $ENV{QUIET}; return $g; }
# irreducible multisets (parts >= 2) with D <= Dmax, as sorted arrays
sub irr_multisets { my ($Dmax) = @_; my @out;
  my $rec; $rec = sub { my ($minj, $D, $cur) = @_; push @out, [@$cur] if @$cur;
    for my $j ($minj..$Dmax+1) { last if $D + $j - 1 > $Dmax; $rec->($j, $D + $j - 1, [@$cur, $j]); } };
  $rec->(2, 0, []); return @out; }
for my $x (rq(2), rq(1,3), rq(-5,7), rq(3,2), rq(7,4), rq(-2,9)) {
  my $c = 1 + $x; my $tl = -$x / $c; my $sig = $tl * (1 + $tl);
  print "x = $x  (c = $c, Lemma-A tau = $tl, sigma = $sig)\n";
  my $tau = s_zero($M); $tau->[$_] = rq((-1)**($_+1)) / rfact($_) for 1..$M;
  my $u = s_scale($tau, 4 * $x / $c**2);
  my $mu = s_scale($u, -1);
  my $sq = s_pow_rat($mu, rq(1,2));  my $q4 = s_pow_rat($mu, rq(1,4));
  my $sstar = s_scale(s_inv(s_add(s_const(1, $M), $sq)), 2 / $c);
  my $dphi = s_int(s_sub(s_const(1, $M), s_scale($sstar, $x)));
  my $zeta = s_mul($tau, s_exp(s_scale($dphi, 1 / $c)));
  my $zt = s_shift_down($zeta);                       # valid to order M-1
  my $rho = s_mul(s_shift_down($tau), s_scale(s_inv(s_mul($q4, s_add(s_const(1, $M), $sq))), rq(2)));
  my @taupow = (s_const(1, $M)); push @taupow, s_mul($taupow[-1], $tau) for 1..$M;
  my (%A, %C);
  for my $j (1..$N+2) {
    my $P = rq(1); $P *= ($j - $_ * $c)**2 for 1..$j-1;            # P_{j-1} = prod_{i=1}^{j-1} (j - i c)^2
    $C{$j} = $P / (rfact($j) * rfact($j-1)) * $c**(4 - 2*$j);
    my $F = s_zero($M); my $a = $j * $x / $c;
    for my $n (0..$M) { $F = s_add($F, s_scale($taupow[$n], poch($a, $n)**2 / (poch(rq($j), $n) * rfact($n)))); }
    my $core = s_scale(s_mul(s_powint($zt, $j), $F), rq((-1)**$j) * $C{$j});   # A_j / t^(j-1), valid to order M-1
    my $Aj = s_zero($N); for my $k ($j-1..$N) { $Aj->[$k] = $core->[$k - $j + 1]; }
    $A{$j} = $Aj; }
  # (1), (2)
  for my $k (0..$N) { my ($s0, $s1) = (rq(0), rq(0)); for my $j (1..$N+1) { $s0 += $A{$j}[$k]; $s1 += $j * $A{$j}[$k]; }
    ok("sum rule t^$k", $s0, $k == 0 ? -$c**2 : rq(0)); ok("Sum j A_j = -c^2 rho, t^$k", $s1, -$c**2 * $rho->[$k]); }
  # (3) against the Mathematica fits
  my $D = parse_mlist($kt5src, { x => $x });
  my ($rhodat, $dil) = ($D->[0], $D->[1]);
  for my $k (0..5) { ok("rho t^$k vs kt5 data", $rho->[$k], $rhodat->[$k]); }
  for my $k (1..5) { my $row = $dil->[$k][0];                      # [lam_0, lam_1, ..., lam_{k+1}] coefficients of [t^k](rho q_t - q_0)
    ok("lam_0 coefficient zero at t^$k (data)", $row->[0], rq(0));
    for my $j (1..$k+1) { ok("a_$j^($k) closed form vs kt5 data", $A{$j}[$k], $row->[$j]); } }
  # (4) Cor. 11.3
  for my $k (1..$N) { ok("top a_".($k+1)."^($k)", $A{$k+1}[$k], rq((-1)**($k+1)) * $C{$k+1});
    ok("sub-top a_$k^($k)", $A{$k}[$k], rq((-1)**$k) * $k * $C{$k} * (1 - $x)**2 / (2 * $c**2)); }
  # (5) C_j = (c^2/j) eps_j(sigma)
  for my $j (1..$N+1) { my $e = rq(1); $e *= ($j*$j*$sig + $_*($j - $_)) for 1..$j-1; $e /= rfact($j-1)**2;
    ok("C_$j = (c^2/$j) eps_$j(sigma)", $C{$j}, $c**2 / $j * $e); }
  # (6) Lemma 5.2, resummed form
  my @Y; $Y[1] = [ @{$A{1}} ]; $Y[1][0] += $c**2; $Y[$_] = $A{$_} for 2..$N+1;
  my ($n0, $n1, $n2) = (0, 0, 0);
  for my $R (irr_multisets($N - 1)) {
    my $J = 0; $J += $_ for @$R; my $DR = $J - scalar(@$R);
    my %m; $m{$_}++ for @$R;
    my $G = s_pow_rat(s_scale($Y[1], -1 / $c**2), rq(-$J));
    my $aut = rq(1); my $lead = rq(1);
    for my $j (keys %m) { $G = s_mul($G, s_powint($Y[$j], $m{$j})); $aut *= rfact($m{$j}); $lead *= (rq((-1)**$j) * $C{$j})**$m{$j}; }
    $G = s_scale($G, 1 / $aut);
    for my $k (0..$DR-1) { ok("G_(@$R) has no t^$k", $G->[$k], rq(0)); }
    $n0++ if ok("[t^D] G_(@$R) = prod (-1)^j C_j/|Aut|", $G->[$DR], $lead / $aut);
    $n1++ if ok("[t^(D+1)] G_(@$R) = 0 (sub-top cancellation)", $G->[$DR+1], rq(0));
    $n2++ if $DR + 2 <= $N && $G->[$DR+2] != 0;
  }
  my @all = irr_multisets($N - 1); my $cnt2 = grep { my $J = 0; $J += $_ for @$_; $J - scalar(@$_) + 2 <= $N } @all;
  print "   Lemma 5.2: ", scalar(@all), " irreducible R with D <= ", $N-1, ": leading ok $n0, sub-top cancels $n1; deficit-2 coefficient non-zero in $n2 of $cnt2\n";
}
print $bad ? "SOME MISMATCH ($bad)\n" : "ALL OK\n";
