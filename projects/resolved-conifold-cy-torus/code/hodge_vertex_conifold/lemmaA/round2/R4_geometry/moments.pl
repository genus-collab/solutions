# moments.pl h [outfile] : one-leaf K-twisted moments with the REFEREE engine (RefHodge.pm, no TR, no kappa classes):
#   S_mu(a) = Sum_m (-1)^m/m! Int_{Mbar_{h,1+m}} ch^mu psi_1^a prod_{k=1}^m Tt(psi_{1+k}),  Tt = psi^2/(1-psi),
# i.e. S_mu(a) = Int_{Mbar_{h,1}} ch^mu K psi^a  (K = exp(-Sum kappa_k/k), dilaton leaves as explicit points),
# for EVERY ch-monomial mu = ch_1^m1 ch_3^m3 ... (odd l <= lmax, lmax = largest odd <= h) of degree <= 3h-3, a = 0..3h-2.
# Output lines: "mu | a:value ; a:value ..." with mu = "m1,m3,m5,m7[,m9]".
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($h, $out) = @ARGV; $out //= "$FindBin::Bin/refdata/moments_h${h}.txt"; my $t0 = time;
my $maxd = $h == 1 ? 1 : 3*$h - 3;
my $NCH = $RefHodge::NCH; my @CHL = @RefHodge::CHL;
my @use = grep { $CHL[$_] <= ($h < 1 ? 1 : $h) } 0..$NCH-1;   # ch_l with l <= h suffice for lambda classes
# enumerate monomials of degree <= maxd in the allowed ch_l
my @mons = ([ (0) x $NCH ]);
for my $i (@use) { my @new; for my $m (@mons) { my $d = 0; $d += $CHL[$_] * $m->[$_] for 0..$NCH-1; for (my $k = 0; $d + $k * $CHL[$i] <= $maxd; $k++) { my @x = @$m; $x[$i] = $k; push @new, \@x; } } @mons = @new; }
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
open my $o, '>', $out or die;
my $cnt = 0;
for my $mu (sort { my $da = 0; my $db = 0; $da += $CHL[$_] * $a->[$_] for 0..$NCH-1; $db += $CHL[$_] * $b->[$_] for 0..$NCH-1; $da <=> $db } @mons) {
  my $dm = 0; $dm += $CHL[$_] * $mu->[$_] for 0..$NCH-1;
  my $B = 3*$h - 2 - $dm; next if $B < 0;   # total psi degree on Mbar_{h,1+m}: a + Sum(D) = 3h-2+m - dm
  my %S;
  for my $P (@{RefHodge::allparts($B)}) { my $m = scalar @$P; my $ex = 0; $ex += $_ for @$P; my $t = $B - $ex;   # leaf exponent a = t
    my @D = map { $_ + 1 } @$P; my %mult; $mult{$_}++ for @D; my $den = Q(1); $den *= fact($_) for values %mult; my $w = (($m % 2) ? Q(-1) : Q(1)) / $den;
    my $v = RefHodge::I($h, $mu, [$t, @D]); next if $v == 0;
    $S{$t} = ($S{$t} // Q(0)) + $w * $v; }
  my @terms = map { "$_:" . $S{$_} } sort { $a <=> $b } grep { $S{$_} != 0 } keys %S;
  print $o join(',', @$mu), " | ", join(' ; ', @terms), "\n"; $cnt++; }
close $o;
printf "genus %d: %d ch-monomials, moments written to %s [%.1f s]\n", $h, $cnt, $out, time - $t0;
