# logfit.pl : given a series C(lambda) = 1 + Sum_k C_k lambda^{2k} (sigma-poly coefficients), test whether
#   log C = Sum_t x_t Sum_k c_k (e1^{2k} e2^t + e2^{2k} e1^t) lambda^{2k},   c_k = -B_{2k}/(2k (2k)!),  e1 + e2 = 1, e1 e2 = -sigma,
# with universal constants x_t (t = 0..T), i.e. C = prod S(e1 lambda)^{x(e2)} S(e2 lambda)^{x(e1)}, x(e) = Sum_t x_t e^t.
# Usage: perl logfit.pl T "c1_sigma_coeffs" "c2..." ...  where each argument is a comma list of sigma coefficients of C_k.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use SerS;
sub R { Math::BigRat->new(shift) }
my ($T, @args) = @ARGV; my $K = scalar @args;
my $C = SerS::szero(); $C->[0] = [R(1)]; for my $k (1..$K) { $C->[$k] = [ map { R($_) } split /,/, $args[$k-1] ]; }
# log of series
my $L = SerS::szero(); { my $X = SerS::sadd($C, SerS::sscale(SerS::sone(), R(-1))); my $pw = SerS::sone();
  for my $n (1..$K) { $pw = SerS::smul($pw, $X); $L = SerS::sadd($L, SerS::sscale($pw, R((($n % 2) ? 1 : -1)) / $n)); } }
# basis g_{t,k} = e1^{2k} e2^t + e2^{2k} e1^t as sigma-polys: use power sums of (e1,e2): e1^a e2^t + e2^a e1^t
sub mix { my ($a, $t) = @_; # e1^a e2^t + e2^a e1^t = (e1e2)^min * (e1^{a-m} e2^{t-m} + ...)
  my $m = $a < $t ? $a : $t; my $p = SerS::pn(abs($a - $t)); my $pr = [R(1)]; for (1..$m) { $pr = SerS::pmul($pr, [R(0), R(-1)]); }
  my $r = SerS::pmul($pr, $p); return $a == $t ? SerS::pscale($r, R(1)/2) : $r; }   # a == t: e1^a e2^a counted twice -> divide (we want e1^a e2^a + e2^a e1^a = 2 (e1e2)^a)
# linear system: for k = 1..K, sigma-coefficients: [L_k] = c_k Sum_t x_t mix(2k, t)
my @rows; my @rhs;
for my $k (1..$K) { my $ck = SerS::ck($k); my @cols = map { my $g = mix(2*$k, $_); $g = SerS::pscale($g, 2) if 2*$k == $_; SerS::pscale($g, $ck) } 0..$T;
  my $lk = $L->[$k]; my $deg = 0; for (@cols, $lk) { $deg = $#$_ if $#$_ > $deg; }
  for my $sp (0..$deg) { push @rows, [ map { $_->[$sp] // R(0) } @cols ]; push @rhs, $lk->[$sp] // R(0); } }
# Gaussian elimination
my $n = $T + 1; my @M = map { [ @{$rows[$_]}, $rhs[$_] ] } 0..$#rows; my $row = 0; my @piv;
for my $c (0..$n-1) { my $p; for my $i ($row..$#M) { if ($M[$i][$c] != 0) { $p = $i; last } } next unless defined $p; @M[$row, $p] = @M[$p, $row];
  my $pv = $M[$row][$c]; $_ = $_ / $pv for @{$M[$row]}; for my $i (0..$#M) { next if $i == $row || $M[$i][$c] == 0; my $f = $M[$i][$c]; $M[$i][$_] -= $f * $M[$row][$_] for 0..$n; } push @piv, $c; $row++; }
my $ok = 1; for my $i ($row..$#M) { $ok = 0 if $M[$i][$n] != 0; }
print "log-fit with t = 0..$T: ", ($ok ? "CONSISTENT" : "inconsistent"), " (rank $row, equations ", scalar(@M), ")\n";
if ($ok) { for my $k (0..$#piv) { print "   x_$piv[$k] = $M[$k][$n]\n"; } }
for my $k (1..$K) { print "   [l^", 2*$k, "] log C = ", join(' ', map { "($L->[$k][$_])s^$_" } grep { ($L->[$k][$_] // 0) != 0 } 0..$#{$L->[$k]}), "\n"; }
