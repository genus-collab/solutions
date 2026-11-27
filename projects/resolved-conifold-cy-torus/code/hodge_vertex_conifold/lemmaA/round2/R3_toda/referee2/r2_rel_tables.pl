#!/usr/bin/perl
# Second referee: the leading-order relations of P_k (k = 2, 3, 4) evaluated EXACTLY with Route E's tables (genus <= 6),
# by the independent engine R2rel.pm, plus an independent re-coding of the author's Theorem M formula (report s.4), which
# must equal 2 x (engine value for P_2 = (D_1^2 + D_2)/2).
# Usage: perl r2_rel_tables.pl k Gmin Gmax rmax Jmax taus   e.g.  perl r2_rel_tables.pl 2 0 4 2 9 "2/7,-5/3"
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use R2rel; use R2tab;
$| = 1;
my ($kP, $Gmin, $Gmax, $rmax, $Jmax, $taus) = @ARGV; $taus //= "2/7,-5/3";
my @taus = map { Math::BigRat->new($_) } split /,/, $taus;
my $P = R2rel::hirota_Pk($kP); print "# P_$kP = ", R2rel::show_P($P), "\n";
my $Q = R2rel::ring_Q();
my %cache;
sub Xv { my ($g, $U, $tau) = @_; my $key = "$g|" . join(',', @$U) . "|$tau"; return $cache{$key} //= R2tab::xi_at($g, $U, $tau); }
# author's Theorem M (report s.4), re-coded from the statement
sub ff { my ($x, $k) = @_; my $p = Math::BigRat->new(1); $p *= ($x - $_) for 0..$k-1; $p }
sub fct { my $n = shift; my $f = Math::BigRat->new(1); $f *= $_ for 2..$n; $f }
sub thmM { my ($G, $eps, $Rp, $tau) = @_; my $n = 2*$G + $eps; my $sig = $tau * ($tau + 1); my $Pi2 = 2*$tau + 1;
  my @R = sort { $a <=> $b } @$Rp; my $J = 0; $J += $_ for @R; my $D = $J - @R; my @R2 = sort { $a <=> $b } (@R, 2); my $D2 = $D + 1;
  my $X = sub { my ($g, $T) = @_; Xv($g, [sort { $a <=> $b } @$T], $tau) };
  my $s = Math::BigRat->new(0);
  for my $g (0..$G+1) { my $k = $n + 1 - 2*$g; next if $k < 1; my $d = $D2 - 2*$g + 2; my $f = ff($d, $k); next if $f == 0; $s -= 2 * $Pi2 * $sig**$g * $tau**$k * $f * $X->($g, \@R2) / fct($k); }
  if ($eps == 0) { my $Bv; if (@R) { $Bv = $J * ($J + 1) * $X->($G, \@R); } else { $Bv = $G == 0 ? Math::BigRat->new('1/2') : ($G == 1 ? -(1 + 1/$sig)/24 : Math::BigRat->new(0)); } $s += 2 * $sig**($G+1) * $Bv; }
  for my $g (0..$G) { my $k = $n - 2*$g; next if $k < 1; my $Cv;
    if (@R) { my $f = ff($D - 2*$g + 2, $k); next if $f == 0; $Cv = $J * ($J + 1) * $f * $X->($g, \@R); } else { $Cv = ($g == 0 && ($k == 1 || $k == 2)) ? Math::BigRat->new(1) : Math::BigRat->new(0); }
    $s += $sig**($g+1) * $tau**$k * $Cv / fct($k); }
  my $Phi = sub { my ($g, $k, $T) = @_; if (!@$T) { return Math::BigRat->new(0) if $g > 0; return $tau if $k == 1; return $tau*$tau/2 if $k == 2; return Math::BigRat->new(0); }
    my $JT = 0; $JT += $_ for @$T; my $f = ff($JT - @$T - 2*$g + 2, $k); return Math::BigRat->new(0) if $f == 0; $sig**$g * $tau**$k * $JT * $f * $X->($g, $T) / fct($k) };
  my $r = @R;
  for my $mask (0..(1 << $r) - 1) { my @T1 = map { $R[$_] } grep { ($mask >> $_) & 1 } 0..$r-1; my @T2 = map { $R[$_] } grep { !(($mask >> $_) & 1) } 0..$r-1;
    for my $g1 (0..$n) { for my $k1 (1..$n + 2) { my $rest = $n + 2 - 2*$g1 - $k1; last if $rest < 1; my $f1 = $Phi->($g1, $k1, \@T1); next if $f1 == 0;
        for my $g2 (0..$rest) { my $k2 = $rest - 2*$g2; last if $k2 < 1; my $f2 = $Phi->($g2, $k2, \@T2); next if $f2 == 0; $s += $f1 * $f2; } } } }
  $s }
my @Rs = ([]); my $recR; $recR = sub { my ($pre, $min, $J) = @_; for (my $j = $min; $J + $j <= $Jmax; $j++) { my @q = (@$pre, $j); push @Rs, [@q]; $recR->(\@q, $j, $J + $j) if @q < $rmax; } }; $recR->([], 2, 0);
my ($n, $bad, $badM, $nM) = (0, 0, 0, 0);
for my $G ($Gmin..$Gmax) { for my $Rp (@Rs) { for my $eps (0, 1) {
  my @res; my $ok = 1; my $okM = 1; my $nt = 0;
  for my $tau (@taus) { my $X = sub { my ($g, $U) = @_; Xv($g, $U, $tau) };
    my ($v, $nterms) = eval { R2rel::relation(P => $P, G => $G, eps => $eps, Rparts => $Rp, ring => $Q, tau => $tau, X => $X) };
    if ($@) { $ok = -1; last; } $nt = $nterms; push @res, $v; $ok = 0 if $v != 0;
    if ($kP == 2) { my $m = eval { thmM($G, $eps, $Rp, $tau) }; if (!$@) { $nM++; if ($m != 2 * $v) { $okM = 0; } } } }
  next if $ok == -1;
  $n++; $bad++ unless $ok; $badM++ unless $okM;
  my $J = 0; $J += $_ for @$Rp; my $Dk = $J - @$Rp + $kP - 1;
  printf "P_%d G=%d eps=%d R=(%s) [R+%d: D=%d excess=%d] terms=%d: %s%s\n", $kP, $G, $eps, join(',', @$Rp), $kP, $Dk, 2*$G - 2 - $Dk, $nt, ($ok ? "0" : "NONZERO " . join(' ', @res)), ($kP == 2 ? ($okM ? "; thmM = 2*engine" : "; thmM DIFFERS") : "");
} } }
printf "SUMMARY P_%d G=%d..%d rmax=%d Jmax=%d: %d relations, %d nonzero%s\n", $kP, $Gmin, $Gmax, $rmax, $Jmax, $n, $bad, ($kP == 2 ? "; author's Theorem M formula vs engine: $badM disagreements in $nM evaluations" : "");
