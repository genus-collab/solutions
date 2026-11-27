#!/usr/bin/perl
# xy_twopoint.pl tau g jmax : two-leaf values from the ABDKS x-y swap formula (LogTR, y affine; arXiv:2312.16950 Prop. 3.1, n = 2),
# after integration by parts, compared with the Hodge-side engine tables of Route E (engine_h<g>_r2.txt).
# For integer tau >= 1, mu = (m1, m2), with hb = hbar, e2 = 1 + tau, Q_mu(z) = prod_{M2(mu)} (z + 1 + m' hb/e2) / prod_{M1(mu)} (z + m hb/tau),
#   [x^mu] F_{g,2} = (1/(m1 m2)) [hb^{2g}] Res_{z2=oo} Res_{z1=oo} Q_m1(z1) Q_m2(z2) dz1 dz2 / ((z1-z2)^2 - k^2 hb^2),  k = (m1+m2)/2,
# iterated residue in the region |z1| > |z2|.  With t = 1/z:  Q_mu = t^{-mu} q_mu(t),
#   iterated residue = Sum_m (k hb)^{2m} Sum_{l + a = m1 - 1 - 2m} C(2m+1+l, l) [t^a] q_m1 * [t^{m2+1+l}] q_m2.
# Dictionary P2: Xi_g(m1, m2) = sigma^{g-1} [x^mu] F_{g,2} / (Pi_m1 Pi_m2).
use strict; use warnings;
use Math::BigRat;
sub Q { Math::BigRat->new(@_) }
my ($tau, $g, $jmax, $datafile) = @ARGV; $tau //= 1; $g //= 2; $jmax //= 6;
$datafile //= "$0"; $datafile =~ s{[^/\\]*$}{}; $datafile .= "../../E_formulas/engine_h${g}_r2.txt";
$tau = Q($tau); my $e2 = $tau + 1; my $sig = $tau * $e2;
my $H = 2*$g;   # hbar truncation
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; return $r; }
sub binomq { my ($n, $k) = @_; return Q(0) if $k < 0 || $k > $n; return fact($n) / (fact($k) * fact($n - $k)); }
# q_mu(t) = prod_{M2} (1 + (1 + m' hb/e2) t) / prod_{M1} (1 + m hb t / tau), as array over t-degree (0..T) of arrays over hb-degree (0..H)
sub qser { my ($mu, $T) = @_; my $n1 = $tau * $mu; my $n2 = $e2 * $mu;
  my @P = map { [ map { Q(0) } 0..$H ] } 0..$T; $P[0][0] = Q(1);
  my $mulf = sub { my ($F) = @_; my @C = map { [ map { Q(0) } 0..$H ] } 0..$T;
    for my $i (0..$T) { for my $k (0..$T-$i) { for my $p (0..$H) { next if $P[$i][$p]->is_zero; for my $q (0..$H-$p) { next if $F->[$k][$q]->is_zero; $C[$i+$k][$p+$q] += $P[$i][$p] * $F->[$k][$q]; } } } }
    @P = @C; };
  for my $kk (1..$n2) { my $m = Q($kk) - ($n2 + 1)/2; my @F = map { [ map { Q(0) } 0..$H ] } 0..$T; $F[0][0] = Q(1); $F[1][0] = Q(1); $F[1][1] = $m / $e2; $mulf->(\@F) if $T >= 1; }
  for my $kk (1..$n1) { my $m = Q($kk) - ($n1 + 1)/2; my @F = map { [ map { Q(0) } 0..$H ] } 0..$T; my $cf = Q(1);
    for my $d (0..$T) { $F[$d][$d] = $cf if $d <= $H; $cf *= -$m / $tau; } $mulf->(\@F); }
  return \@P; }
sub Pi { my $j = shift; my $p = Q(1); $p *= ($j * $tau + $_) for 1..$j-1; return $p / fact($j-1); }
# engine data
my %C; open my $fh, '<', $datafile or die "no $datafile";
while (my $line = <$fh>) { chomp $line; my @f = split /\s*\|\s*/, $line; next unless @f >= 4; my ($h, $r, $k, $rest) = @f; next unless $h == $g && $r == 2;
  for my $term (split /\s*;\s*/, $rest) { next unless $term =~ /(\d+),(\d+):(\S+)/; $C{$k}{"$1,$2"} = Q($3); } }
close $fh;
sub xi_engine { my ($j1, $j2) = @_; my $v = Q(0); my $vv = -Q(1) / $sig;
  for my $k (keys %C) { my $s = Q(0); for my $key (keys %{$C{$k}}) { my ($a, $b) = split /,/, $key; $s += $C{$k}{$key} * Q($j1)**$a * Q($j2)**$b; } $v += $s * $vv**$k; } return $v; }
my ($nb, $nt) = (0, 0);
for my $m1 (1..$jmax) { for my $m2 (1..$m1) {
  my $T = $m2 + 1 + $m1;   # enough t-degree
  my $q1 = qser($m1, $T); my $q2 = qser($m2, $T);
  my $kk = Q($m1 + $m2) / 2;
  my @res = map { Q(0) } 0..$H;   # coefficients in hb
  for (my $m = 0; 2*$m <= $H && 2*$m <= $m1 - 1; $m++) {
    for my $l (0..($m1 - 1 - 2*$m)) { my $a = $m1 - 1 - 2*$m - $l; my $b = $m2 + 1 + $l; next if $b > $T;
      my $bc = binomq(2*$m + 1 + $l, $l);
      for my $p (0..$H) { next if $q1->[$a][$p]->is_zero; for my $q (0..$H - $p - 2*$m) { next if $q2->[$b][$q]->is_zero;
          $res[$p + $q + 2*$m] += $bc * $kk**(2*$m) * $q1->[$a][$p] * $q2->[$b][$q]; } } } }
  my $F = $res[2*$g] / ($m1 * $m2);
  my $xi_tr = $sig**($g-1) * $F / (Pi($m1) * Pi($m2));
  my $xi_h = xi_engine($m1, $m2);
  $nt++; my $ok = ($xi_tr == $xi_h); $nb++ unless $ok;
  printf "tau=%s g=%d (%d,%d): TR-xy %s  engine %s  %s%s\n", $tau, $g, $m1, $m2, $xi_tr, $xi_h, $ok ? "OK" : "MISMATCH", ($m1 + $m2 <= 2*$g - 1 ? "  [range]" : "");
} }
print "total $nt, mismatches $nb\n";
