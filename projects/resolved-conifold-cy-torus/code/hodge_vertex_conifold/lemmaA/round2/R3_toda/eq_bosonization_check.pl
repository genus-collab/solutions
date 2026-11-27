#!/usr/bin/perl
# Check of the explicit form of referee B's operator Eq (line (2) of the task):
#   branching rule at the principal specialisation: d/dp'_1 s_nu(p')|_{q^rho} = eps(nu) s_nu(q^rho),
#   eps(nu) = (q^{1/2}-q^{-1/2}) Sum_{box} q^{-c(box)},
#   and, as an operator on polynomials in p (bosonised -E_0(-x) - 1/vs(x)):
#   Eq f = (1/vs_1) [w^0] ( exp(-Sum_k vs_k w^k p_k / k) f(p_k - vs_k w^{-k}) - f ),   vs_k = q^{k/2} - q^{-k/2}.
# Test: Eq s_nu = eps(nu) s_nu for all |nu| <= 4, at q = Q^2 with rational Q (exact), and the branching rule itself.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use mvlib ();
sub R { Math::BigRat->new(@_) }
my $Q = R($ARGV[0] // "3/2"); my $NMAX = $ARGV[1] // 4;
sub vs { my $k = shift; return $Q ** $k - $Q ** (-$k); }       # q^{k/2} - q^{-k/2}, q = Q^2
# polynomials in p_1..p_N: hash "e1,e2,...,eN" -> BigRat ; w-Laurent: hash wpow -> poly
my $N = $NMAX;
sub key { join(',', @_) } sub unkey { split /,/, $_[0] }
sub padd { my ($a, $b, $c) = @_; $c //= R(1); my %r = %$a; for my $k (keys %$b) { $r{$k} = ($r{$k} // R(0)) + $c * $b->{$k}; delete $r{$k} if $r{$k} == 0; } \%r }
sub pmul { my ($a, $b) = @_; my %r; for my $k1 (keys %$a) { my @e1 = unkey($k1); for my $k2 (keys %$b) { my @e2 = unkey($k2);
      my $k = key(map { $e1[$_] + $e2[$_] } 0..$N-1); $r{$k} = ($r{$k} // R(0)) + $a->{$k1} * $b->{$k2}; } } for (keys %r) { delete $r{$_} if $r{$_} == 0 } \%r }
# Schur polynomial s_nu in power sums: Sum_mu chi^nu(mu) p_mu / z_mu
sub schur { my ($nu) = @_; my $n = 0; $n += $_ for @$nu; my %s;
  for my $mu (mvlib::partitions($n)) { my $ch = mvlib::chi($nu, $mu); next unless $ch; my @e = (0) x $N; $e[$_-1]++ for @$mu;
    $s{key(@e)} = R($ch) / mvlib::zee(@$mu); } \%s }
# Eq f via the bosonic formula
sub Eq { my ($f) = @_;
  # shifted f: substitute p_k -> p_k - vs_k w^{-k}: result as hash wpow -> poly
  my %sh = (0 => {});
  for my $k (keys %$f) { my @e = unkey($k); my %term = (0 => { key((0) x $N) => $f->{$k} });
    for my $i (0..$N-1) { for (1..$e[$i]) {   # multiply by (p_{i+1} - vs w^{-(i+1)})
        my %new; for my $wp (keys %term) { my $pl = $term{$wp};
          my @u = (0) x $N; $u[$i] = 1; my $m1 = pmul($pl, { key(@u) => R(1) }); $new{$wp} = padd($new{$wp} // {}, $m1);
          my $wq = $wp - ($i + 1); $new{$wq} = padd($new{$wq} // {}, $pl, -vs($i + 1)); }
        %term = %new; } }
    for my $wp (keys %term) { $sh{$wp} = padd($sh{$wp} // {}, $term{$wp}); } }
  # multiply by exp(-Sum_k vs_k w^k p_k / k) and take [w^0]: need w^m with m = -wp (wp <= 0)
  my $res = {};
  for my $wp (keys %sh) { my $m = -$wp; next if $m < 0;
    # [w^m] exp(-Sum_k vs_k w^k p_k/k) = Sum over partitions lambda of m: prod (-vs_l p_l / l)^{m_l} / m_l!
    my $coef = {};
    for my $lam (mvlib::partitions($m)) { my %ml; $ml{$_}++ for @$lam; my $c = R(1); my @e = (0) x $N; my $ok = 1;
      for my $l (keys %ml) { if ($l > $N) { $ok = 0; last } $c *= (-vs($l) / $l) ** $ml{$l} / mvlib::fact($ml{$l}); $e[$l-1] += $ml{$l}; }
      next unless $ok; $coef = padd($coef, { key(@e) => $c }); }
    $res = padd($res, pmul($coef, $sh{$wp})); }
  $res = padd($res, $f, R(-1));
  my %out = map { $_ => $res->{$_} / vs(1) } keys %$res; return \%out; }
sub contents { my ($nu) = @_; my @c; for my $i (0..$#$nu) { for my $j (0..$nu->[$i]-1) { push @c, $j - $i; } } @c }
my ($n, $bad) = (0, 0);
for my $sz (1..$NMAX) { for my $nu (mvlib::partitions($sz)) { my $s = schur($nu); my $lhs = Eq($s);
    my $eps = R(0); $eps += $Q ** (-2*$_) for contents($nu); $eps *= vs(1);
    my $diff = padd($lhs, $s, -$eps); my $ok = !%$diff; $n++; $bad++ unless $ok;
    printf "nu=(%s): Eq s_nu = eps s_nu %s\n", join(',', @$nu), ($ok ? "OK" : "FAIL"); } }
# branching rule itself: d/dp'_1 s_nu at p'_k = 1/vs_k equals eps(nu) s_nu(q^rho)
my ($nb, $bb) = (0, 0);
for my $sz (1..$NMAX) { for my $nu (mvlib::partitions($sz)) { my $s = schur($nu); my ($v, $dv) = (R(0), R(0));
    for my $k (keys %$s) { my @e = unkey($k); my $t = $s->{$k}; for my $i (0..$N-1) { $t *= (R(1)/vs($i+1)) ** $e[$i] if $e[$i]; } $v += $t;
      if ($e[0] > 0) { my $t2 = $s->{$k} * $e[0]; for my $i (0..$N-1) { my $ee = $e[$i] - ($i == 0 ? 1 : 0); $t2 *= (R(1)/vs($i+1)) ** $ee if $ee; } $dv += $t2; } }
    my $eps = R(0); $eps += $Q ** (-2*$_) for contents($nu); $eps *= vs(1); $nb++; $bb++ if $dv != $eps * $v; } }
printf "SUMMARY (Q=%s, |nu| <= %d): bosonised Eq: %d Schur functions, %d failures; branching rule: %d, %d failures\n", $Q, $NMAX, $n, $bad, $nb, $bb;
