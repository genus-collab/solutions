#!/usr/bin/perl
# xy_twopoint_diff.pl tau g jmax : two-point values via the difference form of Theorem 1 (n = 2):
#   [x^{m1} x^{m2}] F_{g,2} = (1/(m1 m2)) [hb^{2g}] (1/(2 k hb)) [v^{-1}] Q_{m2}(v) ( Q_{m1,+}(v + k hb) - Q_{m1,+}(v - k hb) ),  k = (m1+m2)/2,
# Q_+ = polynomial part of the Laurent expansion at v = oo.  Compared with Route E's engine_h<g>_r2.txt through P2.
use strict; use warnings; use FindBin; use lib "$FindBin::Bin/lib"; use Math::BigRat; use MP;
my ($tau, $g, $jmax) = @ARGV; $tau = Q($tau // 1); $g //= 2; $jmax //= 5; my $e2 = $tau + 1; my $sig = $tau * $e2;
my $H = 2*$g + 1;
sub fact { my $k = shift; my $r = Q(1); $r *= $_ for 2..$k; return $r; }
sub Pi { my $j = shift; my $p = Q(1); $p *= ($j * $tau + $_) for 1..$j-1; return $p / fact($j-1); }
# Laurent coefficients of Q_mu at oo: Q_mu = v^mu sum_a q_a(hb) v^{-a}
sub qser { my ($mu, $T) = @_; my $n1 = $tau * $mu; my $n2 = $e2 * $mu; my @P = map { [ map { Q(0) } 0..$H ] } 0..$T; $P[0][0] = Q(1);
  my $mulf = sub { my ($F) = @_; my @C = map { [ map { Q(0) } 0..$H ] } 0..$T;
    for my $i (0..$T) { for my $k (0..$T-$i) { for my $p (0..$H) { next if $P[$i][$p]->is_zero; for my $q (0..$H-$p) { next if $F->[$k][$q]->is_zero; $C[$i+$k][$p+$q] += $P[$i][$p] * $F->[$k][$q]; } } } } @P = @C; };
  for my $kk (1..$n2) { my $m = Q($kk) - ($n2 + 1)/2; my @F = map { [ map { Q(0) } 0..$H ] } 0..$T; $F[0][0] = Q(1); if ($T >= 1) { $F[1][0] = Q(1); $F[1][1] = $m / $e2; } $mulf->(\@F); }
  for my $kk (1..$n1) { my $m = Q($kk) - ($n1 + 1)/2; my @F = map { [ map { Q(0) } 0..$H ] } 0..$T; my $cf = Q(1); for my $d (0..$T) { $F[$d][$d] = $cf if $d <= $H; $cf *= -$m / $tau; } $mulf->(\@F); }
  return \@P; }
my %C; open my $fh, '<', "$FindBin::Bin/../../E_formulas/engine_h${g}_r2.txt" or die; while (my $l = <$fh>) { chomp $l; my @f = split /\s*\|\s*/, $l; next unless @f >= 4 && $f[0] == $g;
  for my $t (split /\s*;\s*/, $f[3]) { next unless $t =~ /(\d+),(\d+):(\S+)/; $C{$f[2]}{"$1,$2"} = Q($3); } } close $fh;
sub xe { my ($a1, $a2) = @_; my $v = Q(0); for my $k (keys %C) { my $s = Q(0); for my $key (keys %{$C{$k}}) { my ($a, $b) = split /,/, $key; $s += $C{$k}{$key} * Q($a1)**$a * Q($a2)**$b; } $v += $s * (-Q(1)/$sig)**$k; } $v }
my ($nt, $nb) = (0, 0);
for my $m1 (1..$jmax) { for my $m2 (1..$m1) { my $k = Q($m1 + $m2)/2; my $T = $m1 + $m2 + 2;
  my $q1 = qser($m1, $T); my $q2 = qser($m2, $T);
  # polynomial part of Q_m1: sum_{a<=m1} q1_a v^{m1-a}; evaluate D(v) = Q+(v + k hb) - Q+(v - k hb) as poly in (v, hb): keys "e_v,e_h"
  my %D;
  for my $a (0..$m1) { my $deg = $m1 - $a; for my $p (0..$H) { my $c = $q1->[$a][$p]; next if $c->is_zero;
      # (v + k hb)^deg - (v - k hb)^deg = sum_i C(deg,i) v^{deg-i} (k hb)^i (1 - (-1)^i)
      for my $i (1..$deg) { next unless $i % 2; my $bin = fact($deg) / (fact($i) * fact($deg - $i)); my $key = ($deg - $i) . "," . ($p + $i);
        next if $p + $i > $H; $D{$key} = ($D{$key} // Q(0)) + 2 * $c * $bin * $k**$i; } } }
  # [v^{-1}] Q_m2(v) * v^e : coefficient q2_{m2 + e + 1}
  my @r = map { Q(0) } 0..$H;
  for my $key (keys %D) { my ($e, $p) = split /,/, $key; my $a = $m2 + $e + 1; next if $a > $T; for my $q (0..$H - $p) { $r[$p + $q] += $D{$key} * $q2->[$a][$q]; } }
  # divide by 2 k hb: coefficient of hb^{2g} in result = r[2g+1] / (2k)
  my $F = $r[2*$g + 1] / (2 * $k) / ($m1 * $m2);
  my $xi = $sig**($g-1) * $F / (Pi($m1) * Pi($m2)); my $x2 = xe($m1, $m2);
  $nt++; $nb++ if $xi != $x2; } }
print "tau=$tau g=$g (difference form): total $nt, mismatches $nb\n";
