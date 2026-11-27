#!/usr/bin/perl
# Hurwitz limit (tau -> infinity, top power tau^{n+2}) of the leading-order mKP relation M_{G,eps}(R), n = 2G + eps:
#  -4 Sum_{2g+k-1=n} (d_g(R+2))_k H_g(R+2)/k!  + [eps=0] 2 J(J+1) H_G(R)  + Sum_{2g+k=n, k>=1} J(J+1) (d_g(R))_k H_g(R)/k!
#  + Sum_{S subset R} Sum_{Sum(2g_i+k_i-1)=n} phi(g1,k1,R_S) phi(g2,k2,R_S^c) = 0,
#  phi(g,k,T) = J_T (d_g(T))_k H_g(T)/k!  (T nonempty),  phi(0,1,empty) = 1, phi(0,2,empty) = 1/2,
#  H_g = Xi^(0)_g (Hurwitz component), computed by Route B's recursion (REC) (Rec.pm; valid for delta != 0, all genera).
#  (top coefficients: Pi_2 ~ 2 tau, sigma^g ~ tau^{2g}, Xi_g -> Xi^(0)_g).
# Usage: perl mkp_hurwitz_limit.pl Gmin Gmax rmax  (targets R with D(R) >= 2G-1, so that every needed value has delta != 0)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Rec;
$| = 1;
sub R { Math::BigRat->new(@_) }
my ($Gmin, $Gmax, $rmax, $extra) = @ARGV; $rmax //= 2; $extra //= 3;
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub ff { my ($x, $k) = @_; my $p = R(1); for my $i (0..$k-1) { $p *= ($x - $i); } $p }
sub JD { my $T = shift; my $J = 0; $J += $_ for @$T; ($J, $J - scalar(@$T)) }
sub H { my ($g, $T) = @_; my $v = Rec::xi($g, $T); die "H undefined (delta = 0) for g=$g T=@$T" unless defined $v; $v }
sub phi { my ($g, $k, $T) = @_;
  if (!@$T) { return R(0) if $g > 0; return R(1) if $k == 1; return R("1/2") if $k == 2; return R(0); }
  my ($J, $D) = JD($T); my $f = ff(R($D - 2*$g + 2), $k); return R(0) if $f == 0; return $J * $f * H($g, $T) / fact($k); }
sub subsets { my @T = @_; my @out; my $r = scalar @T;
  for my $mask (0..(1<<$r)-1) { my (@a, @b); for my $i (0..$r-1) { if ($mask & (1<<$i)) { push @a, $T[$i] } else { push @b, $T[$i] } } push @out, [\@a, \@b]; } @out }
sub rel { my ($G, $Rr, $par) = @_; my @Rr = @$Rr; my ($J, $D) = JD(\@Rr); my $n = 2*$G + $par; my @R2 = (@Rr, 2); my ($J2, $D2) = JD(\@R2);
  my $t = R(0);
  for my $g (0..$G+1) { my $k = $n + 1 - 2*$g; next if $k < 1; my $f = ff(R($D2 - 2*$g + 2), $k); next if $f == 0; $t += -4 * $f * H($g, \@R2) / fact($k); }
  $t += 2 * $J * ($J+1) * H($G, \@Rr) if $par == 0;
  for my $g (0..$G) { my $k = $n - 2*$g; next if $k < 1; my $f = ff(R($D - 2*$g + 2), $k); next if $f == 0; $t += $J * ($J+1) * $f * H($g, \@Rr) / fact($k); }
  for my $s (subsets(@Rr)) { my ($T1, $T2) = @$s;
    for my $g1 (0..int(($n+1)/2)) { for my $k1 (1..$n+1-2*$g1) { my $rest = $n + 2 - 2*$g1 - $k1; my $f1 = phi($g1, $k1, $T1); next if $f1 == 0;
        for my $g2 (0..int(($rest-1)/2)) { my $k2 = $rest - 2*$g2; next if $k2 < 1; my $f2 = phi($g2, $k2, $T2); next if $f2 == 0; $t += $f1 * $f2; } } } }
  return $t; }
my ($n, $bad) = (0, 0);
for my $G ($Gmin..$Gmax) {
  my @Rs; my $rec; $rec = sub { my ($pre, $min, $D) = @_; for (my $j = $min; $D + $j - 1 <= 2*$G - 1 + $extra; $j++) { my @p = (@$pre, $j); my $DD = $D + $j - 1;
      push @Rs, [@p] if $DD >= 2*$G - 1; $rec->(\@p, $j, $DD) if scalar(@p) < $rmax; } };
  $rec->([], 2, 0);
  for my $Rr (@Rs) { for my $par (0, 1) { my $v = rel($G, $Rr, $par); $n++; $bad++ if $v != 0;
      my ($J, $D) = JD($Rr); printf "G=%d R={%s} %s [D(R)=%d, D(R+2)=%d]: %s\n", $G, join(',', @$Rr), ($par ? "ODD " : "EVEN"), $D, $D+1, ($v == 0 ? "OK" : "FAIL $v"); } } }
printf "SUMMARY Hurwitz limit: %d relations, %d failures\n", $n, $bad;
