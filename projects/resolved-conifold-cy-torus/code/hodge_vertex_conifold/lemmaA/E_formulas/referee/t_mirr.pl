# t_mirr.pl h [rmax] [rmin] : independent check of M_irr(h) with the referee engine (dilaton leaves as explicit points, no kappa
# classes, no K, no forgetful recursion).  For every multiset R of integers >= 2 with D(R) <= 2h-2 (r between rmin and rmax):
# Xi_h(R; tau) is evaluated exactly at h+2 values of tau with distinct sigma = tau(1+tau); the h+1 components Xi^(k)
# (Xi = Sum_k v^k Xi^(k), v = -1/sigma) are recovered by Lagrange interpolation in v and the extra value is used as a check.
# Lemma A points (D <= 2h-3): all components must vanish.  Control points (D = 2h-2): the non-zero components are listed.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($h, $rmax, $rmin) = @ARGV; $rmax //= 99; $rmin //= 1; my $t0 = time;
my @taus = map { Q($_) } qw(1 2 3 1/2 1/3 -1/3 3/2 5 -7/2 7 1/5 -2/7);
@taus = @taus[0..$h+1];
my $maxd = $h == 1 ? 1 : 3*$h - 3;
my @om = map { RefHodge::omega($h, [Q(1), $_, -1 - $_], $maxd) } @taus;
my @vs = map { Q(-1) / ($_ * (1 + $_)) } @taus;
printf "genus %d: %d tau values, Omega built [%.1f s]\n", $h, scalar(@taus), time - $t0;
my @pts; my $rec; $rec = sub { my ($pre, $min, $D) = @_; push @pts, [@$pre] if @$pre && @$pre >= $rmin && @$pre <= $rmax;
  return if @$pre >= $rmax; for (my $j = $min; ($j - 1) + $D <= 2*$h - 2; $j++) { $rec->([@$pre, $j], $j, $D + $j - 1); } };
$rec->([], 2, 0);
# Lagrange interpolation: components from values at v_0..v_h
sub interp { my ($v, $y, $d) = @_; my @C = map { Q(0) } 0..$d;
  for my $i (0..$d) { my @L = (Q(1)); my $den = Q(1);
    for my $m (0..$d) { next if $m == $i; my @N = map { Q(0) } 0..$#L+1; for my $k (0..$#L) { $N[$k+1] += $L[$k]; $N[$k] -= $L[$k] * $v->[$m]; } @L = @N; $den *= ($v->[$i] - $v->[$m]); }
    for my $k (0..$d) { $C[$k] += $y->[$i] * $L[$k] / $den; } }
  return @C; }
my ($nA, $badA, $nC, $nzC, $badI) = (0, 0, 0, 0, 0);
for my $p (sort { scalar(@$a) <=> scalar(@$b) || join(',', @$a) cmp join(',', @$b) } @pts) {
  my $D = 0; $D += $_ - 1 for @$p;
  my @y = map { RefHodge::xi_value(RefHodge::xi($h, $om[$_], $p), Q(-1)) } 0..$#taus;
  my @C = interp(\@vs, \@y, $h);
  my $chk = Q(0); for my $k (0..$h) { $chk += $C[$k] * $vs[$h+1] ** $k; } my $iok = $chk == $y[$h+1]; $badI++ unless $iok;
  my @nz = grep { $C[$_] != 0 } 0..$h;
  if ($D <= 2*$h - 3) { $nA++; $badA++ if @nz || grep { $_ != 0 } @y;
    printf "genus %d R = {%s} D = %d: %s%s [%.1f s]\n", $h, join(',', @$p), $D, ((@nz || grep { $_ != 0 } @y) ? "NON-ZERO components @nz" : "all components 0"), ($iok ? "" : " (INTERPOLATION CHECK FAILED)"), time - $t0; }
  else { $nC++; $nzC++ if @nz; printf "genus %d R = {%s} D = %d (control): non-zero components {%s}%s [%.1f s]\n", $h, join(',', @$p), $D, join(' ', @nz), ($iok ? "" : " (INTERPOLATION CHECK FAILED)"), time - $t0; } }
printf "SUMMARY genus %d (r = %d..%s): %d Lemma-A points, %d with a non-zero component; %d control points, %d non-zero; %d interpolation-check failures; memo %d; %.1f s\n",
  $h, $rmin, ($rmax == 99 ? "all" : $rmax), $nA, $badA, $nC, $nzC, $badI, scalar(keys %RefHodge::IM), time - $t0;
