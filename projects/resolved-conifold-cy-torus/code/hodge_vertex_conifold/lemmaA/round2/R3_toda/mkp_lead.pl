#!/usr/bin/perl
# Leading order at the dilaton point (w = 1 + s -> 0) of the first modified-KP equation of the 2D Toda hierarchy
#   2 d_{p2}(F - F1) + d_{p1}^2 (F + F1) + (d_{p1}(F - F1))^2 = 0,   F = log Z_MV,  F1(p) = F(e^beta o p),  beta = tau x,
# in the variables q_k = x sigma Pi_k p_k (unit legs resummed in y = q_1), sector e^{u(J+2)}, w-order = (x-power) + r.
# With A = F1 - F = Sum_{k>=1} beta^k E^k F / k!:   -2 d_{p2} A + d_{p1}^2 (2F + A) + (d_{p1} A)^2 = 0.
# Leading coefficients: G^(g)_T ~ e^{uJ} X_g(T) w^{2g-2+r},  E^k -> (delta)_k (falling), d_y -> J,  delta_g(T) = D(T) - 2g + 2.
# Usage: perl mkp_lead.pl Gmin Gmax rmax Jmax [taus]   (prints EVEN/ODD residuals; all must be 0)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use XiData;
sub R { Math::BigRat->new(@_) }
my ($Gmin, $Gmax, $rmax, $Jmax, $taustr) = @ARGV;
$Gmin //= 0; $Gmax //= 3; $rmax //= 2; $Jmax //= 8; $taustr //= "2/7,3/2,-5/3";
my @taus = map { R($_) } split /,/, $taustr;
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub ff { my ($x, $k) = @_; my $p = R(1); for my $i (0..$k-1) { $p *= ($x - $i); } $p }
sub Pi { my ($k, $tau) = @_; my $p = R(1); for my $a (1..$k-1) { $p *= ($k*$tau + $a); } $p / fact($k-1) }
sub JD { my $T = shift; my $J = 0; $J += $_ for @$T; return ($J, $J - scalar(@$T)); }
my $NEEDDATA = 0;
my ($PG, $PT, $PD) = split /\|/, ($ENV{PERTURB} // "");   # sensitivity test: PERTURB="g|t1,t2,..|delta"
sub X { my ($g, $T, $tau) = @_; my $v = XiData::xi_at($g, $T, $tau);
  if (defined $PD && $g == $PG && join(',', sort { $a <=> $b } @$T) eq $PT) { $v += R($PD); }
  return $v; }
# factor of (d_{p1} A): (g, k, T)
sub Phi { my ($g, $k, $T, $tau, $sig) = @_;
  if (!@$T) { return R(0) if $g > 0; return $tau if $k == 1; return $tau*$tau/2 if $k == 2; return R(0); }
  my ($J, $D) = JD($T); my $d = $D - 2*$g + 2; my $f = ff(R($d), $k); return R(0) if $f == 0;
  return ($sig ** $g) * ($tau ** $k) * $J * $f * X($g, $T, $tau) / fact($k); }
sub subsets { my @T = @_; my @out; my $r = scalar @T;
  for my $mask (0..(1<<$r)-1) { my (@a, @b); for my $i (0..$r-1) { if ($mask & (1<<$i)) { push @a, $T[$i] } else { push @b, $T[$i] } } push @out, [\@a, \@b]; }
  return @out; }
sub relation { my ($G, $Rr, $tau, $par) = @_;   # par = 0: x^{2G}; 1: x^{2G+1}
  my $sig = $tau * ($tau + 1); my $P2 = Pi(2, $tau); my @Rr = @$Rr; my ($J, $D) = JD(\@Rr);
  my $n = 2*$G + $par; my %part;
  my @R2 = (@Rr, 2); my ($J2, $D2) = JD(\@R2);
  # (a)  -2 d_{p2} A : x-power 2g + k - 1 = n
  my $a = R(0);
  for my $g (0..$G+1) { my $k = $n + 1 - 2*$g; next if $k < 1;
    my $f = ff(R($D2 - 2*$g + 2), $k); next if $f == 0;
    $a += -2 * $P2 * ($sig ** $g) * ($tau ** $k) * $f * X($g, \@R2, $tau) / fact($k); }
  # (b)  2 d_{p1}^2 F : x-power 2g = n
  my $b = R(0);
  if ($par == 0) { my $g = $G;
    if (@Rr) { $b = 2 * ($sig ** ($g+1)) * $J * ($J+1) * X($g, \@Rr, $tau); }
    else { my $B = $g == 0 ? R("1/2") : ($g == 1 ? -(1 + 1/$sig)/24 : R(0)); $b = 2 * ($sig ** ($g+1)) * $B; } }
  # (c)  d_{p1}^2 A : x-power 2g + k = n, k >= 1
  my $c = R(0);
  for my $g (0..$G) { my $k = $n - 2*$g; next if $k < 1;
    if (@Rr) { my $f = ff(R($D - 2*$g + 2), $k); next if $f == 0;
      $c += ($sig ** ($g+1)) * ($tau ** $k) * $J * ($J+1) * $f * X($g, \@Rr, $tau) / fact($k); }
    else { next unless $g == 0; next unless $k == 1 || $k == 2; $c += ($sig ** ($g+1)) * ($tau ** $k) / fact($k); } }
  # (d)  (d_{p1} A)^2 : Sum (2 g_i + k_i - 1) = n
  my $dd = R(0);
  for my $s (subsets(@Rr)) { my ($T1, $T2) = @$s;
    for my $g1 (0..int(($n+1)/2)) { for my $k1 (1..$n+1-2*$g1) { my $rest = $n + 2 - 2*$g1 - $k1;   # = 2 g2 + k2
        my $f1 = Phi($g1, $k1, $T1, $tau, $sig); next if $f1 == 0;
        for my $g2 (0..int(($rest-1)/2)) { my $k2 = $rest - 2*$g2; next if $k2 < 1;
          my $f2 = Phi($g2, $k2, $T2, $tau, $sig); next if $f2 == 0; $dd += $f1 * $f2; } } } }
  return ($a + $b + $c + $dd, $a, $b, $c, $dd);
}
# multisets of parts >= 2 with r <= rmax, J <= Jmax (plus the empty one)
my @Rs = ([]); my $rec; $rec = sub { my ($pre, $min, $J) = @_; for (my $j = $min; $J + $j <= $Jmax; $j++) { my @p = (@$pre, $j); push @Rs, [@p]; $rec->(\@p, $j, $J + $j) if scalar(@p) < $rmax; } };
$rec->([], 2, 0);
my ($ntot, $nbad) = (0, 0);
for my $G ($Gmin..$Gmax) { for my $Rr (@Rs) { for my $par (0, 1) {
  my @R2 = (@$Rr, 2); my ($J2, $D2) = JD(\@R2); my $exc = 2*$G - 2 - $D2;
  my @res; my $ok = 1;
  my $terms = "";
  for my $tau (@taus) { my ($tot, $a, $b, $c, $d) = relation($G, $Rr, $tau, $par); push @res, $tot; $ok = 0 if $tot != 0;
    $terms = "   (a)=$a (b)=$b (c)=$c (d)=$d at tau=$tau" if $ENV{VERBOSE} && $terms eq ""; }
  $ntot++; $nbad++ unless $ok;
  print "$terms\n" if $ENV{VERBOSE};
  printf "G=%d R={%s} %s  [R+2: D=%d excess=%d]  %s\n", $G, join(',', @$Rr), ($par ? "ODD " : "EVEN"), $D2, $exc, ($ok ? "OK (0)" : "FAIL: " . join(' ', @res));
} } }
printf "SUMMARY: %d relations, %d failures\n", $ntot, $nbad;
