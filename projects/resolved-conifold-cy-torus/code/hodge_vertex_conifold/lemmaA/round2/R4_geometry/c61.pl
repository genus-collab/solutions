# c61.pl h1 h2 ... : formula (6.1) against the component tables refdata/ref_engine_h<h>_r1.txt (referee engine, no TR):
#   sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j - l) [lambda^{2h}] S(e1 lambda)^{1 + j e2} S(e2 lambda)^{1 + j e1},
#   e1 = -tau, e2 = 1 + tau, S(y) = (y/2)/sinh(y/2),  Xi = Sum_k v^k Xi^(k), v = -1/sigma.
# Both sides are polynomials in j (deg <= 3h-2) and v (deg <= h); compared at j = 0..3h+1 and h+3 values of tau.
# Also prints, component by component, the comparison of sigma-coefficients (exact polynomial identity in j).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
my $dir = $ENV{DATADIR} // "$FindBin::Bin/refdata";
for my $h (@ARGV) {
  my @C; open my $fh, '<', "$dir/ref_engine_h${h}_r1.txt" or die "no table for h=$h"; while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4; $rest //= '';
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; $C[$k][$m] = Q($cf); } } close $fh;
  my $K = $h; my @sh = map { Q(1) / (Q(4) ** $_ * fact(2*$_ + 1)) } 0..$K;
  my @lg = (Q(0)); for my $n (1..$K) { my $s = $sh[$n] * $n; for my $m (1..$n-1) { $s -= $m * $lg[$m] * $sh[$n - $m]; } push @lg, $s / $n; }
  my @logS = map { -$_ } @lg;
  my ($n, $bad) = (0, 0);
  for my $tau (map { Q($_) } qw(1/3 2 -3/7 5 -1/4 7/2 1/9 3 -2/9 4/5 11/3 -5/2)[0..$h+2]) { my $e1 = -$tau; my $e2 = 1 + $tau; my $sig = $tau * (1 + $tau); my $v = Q(-1) / $sig;
    for my $j (0..3*$h + 1) { my $jq = Q($j);
      my @c = (Q(0)); for my $k (1..$K) { push @c, (1 + $jq * $e2) * $logS[$k] * $e1 ** (2*$k) + (1 + $jq * $e1) * $logS[$k] * $e2 ** (2*$k); }
      my @ex = (Q(1)); for my $nn (1..$K) { my $s = Q(0); for my $m (1..$nn) { $s += $m * $c[$m] * $ex[$nn - $m]; } $ex[$nn] = $s / $nn; }
      my $rhs = $ex[$h]; $rhs *= ($jq - $_) for 1..(2*$h - 2);
      my $lhs = Q(0); for my $k (0..$h) { my $ck = Q(0); for my $a (0..$#{$C[$k] // []}) { next unless defined $C[$k][$a]; $ck += $C[$k][$a] * $jq ** $a; } $lhs += $ck * $v ** $k; }
      $lhs *= $sig ** $h; $n++; $bad++ if $lhs != $rhs; } }
  printf "genus %d: (6.1) vs table: %d comparisons (%d tau x %d j), %d mismatches\n", $h, $n, $h + 3, 3*$h + 2, $bad; }
