# t_basic.pl : sanity checks of the referee engine (RefWK, RefHodge) against classical closed formulas.
#  WK: <tau_{3g-2}>_g = 1/(24^g g!), <tau_0 tau_1 ... > known small values.
#  lambda_g formula   Int_{Mbar_{g,1}} psi^{2g-2} lambda_g = b_g = [t^{2g}] (t/2)/sin(t/2)
#  Faber              Int_{Mbar_g} lambda_g lambda_{g-1} lambda_{g-2} = |B_{2g-2}| |B_{2g}| / (2 (2g-2)! (2g-2) (2g))
#  Getzler-Pandharipande  Int_{Mbar_{g,n}} psi^a lambda_g lambda_{g-1} = (2g-3+n)! |B_{2g}| / (2^{2g-1} (2g)! prod (2a_i-1)!!)
#  Mumford relation   Int_{Mbar_g} lambda_{g-1}^3 = 2 Int lambda_g lambda_{g-1} lambda_{g-2}  (from c(E)c(E^vee) = 1)
#  Harer-Zagier (via Gauss-Bonnet for T(-log)): Int_{Mbar_{g,n}} Lambda^vee(1) K = chi(M_{g,n}) -- tested in t_hz.pl
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefWK; use RefHodge; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
my $gmax = shift // 4; my $t0 = time;
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; $r }
sub df { my $n = shift; my $r = Q(1); for (my $k = $n; $k > 1; $k -= 2) { $r *= $k } $r }
# WK
for my $g (1..6) { my $v = RefWK::wk($g, 3*$g - 2); my $e = Q(1) / (Q(24) ** $g * fact($g)); print "WK <tau_", 3*$g-2, ">_$g ", ($v == $e ? "OK" : "FAIL $v vs $e"), "\n"; }
print "WK <tau_2 tau_3>_2 = ", RefWK::wk(2, 2, 3), " (expected 29/5760)\n";
print "WK <tau_1 tau_1 tau_1>_1... <tau_2 tau_2 tau_2 tau_2>_3 = ", RefWK::wk(3, 2,2,2,2,0), " ; <tau_4 tau_4>_3 = ", RefWK::wk(3, 4, 4), " (expected 607/1451520)\n";
# b_g from (t/2)/sin(t/2) = exp(Sum_k |B_2k| ... ): compute series directly: sin(t/2)/(t/2) = Sum (-1)^k (t/2)^{2k}/(2k+1)!
my $K = 8; my @sn = map { ((($_ % 2) ? Q(-1) : Q(1)) / (Q(4) ** $_) / fact(2*$_ + 1)) } 0..$K; my @inv = (Q(1));
for my $n (1..$K) { my $s = Q(0); for my $m (1..$n) { $s -= $sn[$m] * $inv[$n - $m]; } $inv[$n] = $s; }   # (t/2)/sin(t/2) coefficients of t^{2n}
for my $g (1..$gmax) { my $h = $g; my $maxd = 3*$g - 3 + 1; $maxd = 1 if $g == 1;
  my @lam = RefHodge::lambdas($g, $maxd < 3*$g ? 3*$g : $maxd);
  my $int = sub { my ($cls, $e) = @_; my $s = Q(0); for my $k (keys %$cls) { my @mu = split /,/, $k; $s += $cls->{$k} * RefHodge::I($g, \@mu, $e); } return $s; };
  my $v3 = $int->($lam[$g], [2*$g - 2]); print "genus $g: lambda_g formula ", ($v3 == $inv[$g] ? "OK" : "FAIL $v3 vs $inv[$g]"), "\n";
  if ($g >= 2) { my $c = RefHodge::cmul(RefHodge::cmul($lam[$g], $lam[$g-1], 3*$g), $lam[$g-2], 3*$g); my $v1 = $int->($c, []);
    my $e1 = abs(RefHodge::bern(2*$g - 2)) * abs(RefHodge::bern(2*$g)) / (2 * fact(2*$g - 2) * (2*$g - 2) * (2*$g));
    print "genus $g: Faber lambda_g lambda_{g-1} lambda_{g-2} ", ($v1 == $e1 ? "OK" : "FAIL $v1 vs $e1"), "\n";
    my $c3 = RefHodge::cmul(RefHodge::cmul($lam[$g-1], $lam[$g-1], 3*$g), $lam[$g-1], 3*$g); my $v4 = $int->($c3, []);
    print "genus $g: lambda_{g-1}^3 = 2 x Faber ", ($v4 == 2 * $e1 ? "OK" : "FAIL $v4 vs " . 2*$e1), "\n";
    my $c2 = RefHodge::cmul($lam[$g], $lam[$g-1], 3*$g);
    for my $a ([$g - 1], [$g - 2, 1], [1, 1, $g - 2], [2, $g - 1, 0]) { my $n = scalar @$a; my $sa = 0; $sa += $_ for @$a; next if $sa != $g - 2 + $n;
      my $v2 = $int->($c2, $a); my $e2 = fact(2*$g - 3 + $n) * abs(RefHodge::bern(2*$g)) / (Q(2) ** (2*$g - 1) * fact(2*$g)); $e2 /= df(2*$_ - 1) for @$a;
      print "genus $g: GP psi^{@$a} lambda_g lambda_{g-1} ", ($v2 == $e2 ? "OK" : "FAIL $v2 vs $e2"), "\n"; } } }
printf "time %.1f s, memo %d\n", time - $t0, scalar(keys %RefHodge::IM);
