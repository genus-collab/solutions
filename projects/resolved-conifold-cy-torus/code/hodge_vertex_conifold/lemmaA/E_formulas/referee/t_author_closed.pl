# t_author_closed.pl h : closed formulas NOT used by the author, applied to the AUTHOR's engine (../Hodge.pm, ../HodgePt.pm):
#  (HZ)  Harer-Zagier + Gauss-Bonnet for the log tangent bundle: Xi^(0)_h(0^r) = Int_{Mbar_{h,r}} c(T(-log)) = chi(M_{h,r})
#        = B_{2h}/(2h(2h-2)) prod_{i=0}^{r-1} (2-2h-i)   (uses (L) Lambda^vee(1) K = c(T(-log)) and K; tests the kappa machinery)
#  (G)   lambda_g component Xi^(h)_h(R) = b_h (J-1)...(J-N), N = 2h-3+r, at points inside AND outside the Lemma-A range
#  (MV)  one-part Marino-Vafa formula for Int_{Mbar_{h,1}} Omega lam_d (no K): tests omega_ch and HK with all v-powers.
# Run with HODGE_NCH=4 for genus 7.
use strict; use warnings; use FindBin; use lib "$FindBin::Bin/.."; use Math::BigRat; use Hodge; use HodgePt; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my $h = shift; my $t0 = time;
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; $r }
my @B = (Q(1)); for my $m (1..2*$h + 2) { my $s = Q(0); for my $k (0..$m-1) { my $b = Q(1); for my $i (1..$k) { $b = $b * ($m + 1 - $k + $i) / $i } $s += $b * $B[$k]; } push @B, -$s / ($m + 1); }
Hodge::setup($h);
# (HZ)
my $chi = $B[2*$h] / (2*$h * (2*$h - 2));
{ my $t = Hodge::T(); printf "(HZ) genus %d: Int_{Mbar_h} Lambda^vee(1) K = %s, chi(M_h) = %s : %s [%.1f s]\n", $h, $t->[0], $chi, ($t->[0] == $chi ? "OK" : "FAIL"), time - $t0;
  for my $r (1..3) { my $v = HodgePt::WP([], [(0) x $r]); my $e = $chi; $e *= (2 - 2*$h - $_) for 0..$r-1;
    printf "(HZ) genus %d r = %d: Xi^(0)(0^r) = chi(M_{h,r}) : %s [%.1f s]\n", $h, $r, ((($v->[0] // 0) == $e) ? "OK" : "FAIL $v->[0] vs $e"), time - $t0; } }
# (G)
my $K = $h + 1; my @sn = map { ((($_ % 2) ? Q(-1) : Q(1)) / (Q(4) ** $_) / fact(2*$_ + 1)) } 0..$K; my @inv = (Q(1));
for my $n (1..$K) { my $s = Q(0); for my $m (1..$n) { $s -= $sn[$m] * $inv[$n - $m]; } $inv[$n] = $s; }
my $bh = $inv[$h];
{ my ($n, $bad) = (0, 0);
  for my $R ([2], [5], [2*$h - 1], [2*$h + 3], [0], [1], [2, 2], [3, 2*$h], [0, 7], [2, 2, 2], [1, 3, 2*$h + 1]) { my $v = HodgePt::WP([], $R); my $J = 0; $J += $_ for @$R; my $N = 2*$h - 3 + scalar(@$R);
    my $e = $bh; $e *= ($J - $_) for 1..$N; $n++; $bad++ if ($v->[$h] // 0) != $e; }
  printf "(G) genus %d: lambda_g component = b_h (J-1)...(J-N) at %d points (in and out of range): %d failures [%.1f s]\n", $h, $n, $bad, time - $t0; }
# (MV) without K
{ my $D = 3*$h - 2; my $OM = Hodge::omega_ch($h, $D); my ($n, $bad) = (0, 0);
  my $KK = $h + 1;
  my $fser = sub { my $x = shift; [ map { ((($_ % 2) ? Q(-1) : Q(1)) * $x ** (2*$_)) / (Q(4) ** $_ * fact(2*$_ + 1)) } 0..$KK ] };
  my $smul = sub { my ($a, $b) = @_; my @c = map { Q(0) } 0..$KK; for my $i (0..$KK) { for my $j (0..$KK-$i) { $c[$i+$j] += $a->[$i] * $b->[$j]; } } \@c };
  my $sinv = sub { my $a = shift; my @r = (Q(1) / $a->[0]); for my $nn (1..$KK) { my $s = Q(0); for my $m (1..$nn) { $s += $a->[$m] * $r[$nn-$m]; } push @r, -$s / $a->[0]; } \@r };
  my @hk; for my $a (0..$D) { my $s = [Q(0)]; for my $mk (keys %$OM) { my @m = split /,/, $mk; my $v = Hodge::HK($h, \@m, [], [$a]); next if $v == 0; $s = Hodge::vadd($s, Hodge::vscale($OM->{$mk}, $v)); } push @hk, $s; }
  for my $tau (map { Q($_) } qw(2 -1/3 3/7)) { my $sig = $tau * (1 + $tau); my $vv = Q(-1) / $sig;
    for my $d (1..6) { my $val = Q(0); for my $a (0..$D) { my $c = Q(0); my $p = $hk[$a]; for my $k (0..$#$p) { $c += $p->[$k] * $vv ** $k; } $val += $c * Q($d) ** $a; }
      my $num = [Q(1), map { Q(0) } 1..$KK]; my $den = [Q(1), map { Q(0) } 1..$KK];
      for my $aa (1..$d-1) { $num = $smul->($num, $fser->($d * $tau + $aa)); } for my $aa (1..$d) { $den = $smul->($den, $fser->(Q($aa))); }
      my $pred = $smul->($num, $sinv->($den))->[$h] / ($d * $d);
      $n++; $bad++ if $val * ((-$sig) ** $h) != $pred; } }
  printf "(MV) genus %d: author's Omega/HK vs one-part Marino-Vafa formula: %d comparisons, %d failures [%.1f s]\n", $h, $n, $bad, time - $t0; }
