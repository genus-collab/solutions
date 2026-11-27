# excess_split.pl h : split the one-leaf components Xi^(k)_h(j) = Int c(F-P) p_j(psi) e_k(alpha+alpha^2) (F-P = -R pi_* S_D, D = j x)
# into the EXCESS part (Chern degree m of c(F-P) above its rank N-j: localized on the jump locus P(ker phi), Route C) and the
# NON-EXCESS part (m <= N-j), N = 3h-2.  c(F-P) = Lambda^vee(1) K / prod_{i=1}^j (1 - i psi).  Values at integer j.
# Check: E + NE must equal Xi^(k)_h(j) of the component table.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Pair; use RefHodge;
sub Q { Math::BigRat->new(shift) }
my $h = shift // 3; my $jmax = shift // 2*$h; my $N = 3*$h - 2;
my $M = Pair::load($h);
my @lam = Pair::lambdas($h);
my @L1; for my $i (0..$h) { $L1[$i] = Pair::scal($lam[$i], ($i % 2) ? Q(-1) : Q(1)); }          # degree-i part of c(E^vee) = Lambda^vee(1)
my @Ek = Pair::ek_beta($h); my @Hd = map { Pair::bydeg($_) } @Ek;
# sanity: Lambda^vee(1) e_k = Omega^(k)
my @comp = Pair::omega_components($h);
for my $k (0..$h) { my $P = {}; $P = Pair::add($P, Pair::mul($h, $L1[$_], $Ek[$k])) for 0..$h; my $p1 = Pair::jpoly($h, $P); my $p2 = Pair::jpoly($h, $comp[$k]);
  my $ok = 1; for my $a (0..($#$p1 > $#$p2 ? $#$p1 : $#$p2)) { $ok = 0 if ($p1->[$a] // Q(0)) != ($p2->[$a] // Q(0)); } print "k=$k: Lambda^vee(1) e_k vs Omega^(k): ", ($ok ? "equal" : "DIFFERENT"), "\n"; }
sub mom { my ($C, $a) = @_; my $s = Q(0); for my $mu (keys %$C) { my $v = $M->{$mu}; next unless $v && defined $v->[$a]; $s += $C->{$mu} * $v->[$a]; } return $s; }
sub hcomp { my ($d, $j) = @_; my @c = (Q(1)); for my $i (1..$j) { my @n = map { Q(0) } 0..$d; for my $x (0..$d) { next unless defined $c[$x]; my $p = Q(1); for my $y (0..$d - $x) { $n[$x+$y] += $c[$x] * $p; $p *= $i; } } @c = @n; } return $c[$d] // Q(0); }   # h_d(1..j)
sub ecomp { my ($d, $n) = @_; my @c = (Q(1)); for my $i (1..$n) { my @nn = (@c, Q(0)); for my $x (reverse 1..$#nn) { $nn[$x] += $c[$x-1] * $i if defined $c[$x-1]; } @c = @nn; } return $c[$d] // Q(0); }   # e_d(1..n)
printf "genus %d, N = %d.  k, j, EXCESS part E, NON-EXCESS part NE, total\n", $h, $N;
for my $k (0..$h) { my $tot = Pair::jpoly($h, $comp[$k]);
  for my $j (2..$jmax) {
    my ($E, $NE) = (Q(0), Q(0));
    for my $d1 (0..$h) { for my $d2 (keys %{$Hd[$k]}) { my $P = Pair::mul($h, $L1[$d1], $Hd[$k]{$d2}); next unless %$P;
        for my $a1 (0..$N) { my $c1 = hcomp($a1, $j); next if $c1 == 0;
          for my $a2 (0..$j-1) { my $c2 = ecomp($a2, $j - 1) * (($a2 % 2) ? -1 : 1); next if $c2 == 0;
            my $a = $a1 + $a2; next if $a > $N; my $dK = $N - $d1 - $d2 - $a; next if $dK < 0;
            my $v = mom($P, $a); next if $v == 0; my $m = $d1 + $dK + $a1;
            if ($m > $N - $j) { $E += $c1 * $c2 * $v } else { $NE += $c1 * $c2 * $v } } } } }
    my $t = Pair::peval($tot, Q($j));
    printf "k=%d j=%2d  E=%s  NE=%s  tot=%s%s\n", $k, $j, $E, $NE, $E + $NE, ($E + $NE == $t ? '' : "  MISMATCH table=$t"); } }
