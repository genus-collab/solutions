# check_genus3_alltau.pl -- M_irr(3) for every tau, computer-verified (Perl, exact, no Mathematica, no TR):
#   sigma^3 Xi_3(R) is a polynomial of degree <= 3 in sigma = tau(1+tau) (proof_status s.1.1(e), 1.3), so vanishing at
#   4 values of tau with distinct sigma proves Xi_3(R) = 0 identically. Also prints the tau-free components
#   Xi^(k)_3(R) = Int c(T(-log)) e_k(alpha + alpha^2) prod lam  (k = 0..3), obtained by solving the 4x4 Vandermonde system.
use strict; use warnings; use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }
my $t0 = time;
my @taus = (Q(1), Q(2), Q(-1,3), Q(3,4), Q(2,7));
for my $R ([2],[3],[4],[2,2],[2,3],[2,2,2],[5],[3,3],[2,4],[2,2,3]) {
  my $D = 0; $D += $_ - 1 for @$R;
  my @vals = map { xi(3, $_, @$R) } @taus;
  my @sig = map { $_ * (1 + $_) } @taus;
  # components: sigma^3 Xi = sum_k (-1)^k sigma^(3-k) Xi^(k); solve with the first 4 taus, check the 5th
  my @A = map { my $s = $sig[$_]; [ map { Q((-1)**$_) * $s**(3 - $_) } 0..3 ] } 0..3;
  my @b = map { $sig[$_]**3 * $vals[$_] } 0..3;
  # Gaussian elimination
  my $n = 4; for my $i (0..$n-1) { my $p = $i; $p++ while $A[$p][$i] == 0; @A[$i,$p] = @A[$p,$i]; @b[$i,$p] = @b[$p,$i];
    for my $k (0..$n-1) { next if $k == $i; my $f = $A[$k][$i] / $A[$i][$i]; next if $f == 0; $A[$k][$_] -= $f * $A[$i][$_] for 0..$n-1; $b[$k] -= $f * $b[$i]; } }
  my @X = map { $b[$_] / $A[$_][$_] } 0..$n-1;
  my $pred = 0; $pred += Q((-1)**$_) * $sig[4]**(3 - $_) * $X[$_] for 0..3; my $chk = ($pred == $sig[4]**3 * $vals[4]) ? "5th tau ok" : "5th tau MISMATCH";
  printf "R=(%s) D=%d  Xi at tau=1,2,-1/3,3/4,2/7: %s | Xi^(0..3) = %s | %s [%d s]\n", join(",", @$R), $D, join(" ", @vals), join(", ", @X), $chk, time - $t0;
}
