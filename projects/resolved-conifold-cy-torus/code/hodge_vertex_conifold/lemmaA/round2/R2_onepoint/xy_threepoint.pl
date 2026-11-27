#!/usr/bin/perl
# xy_threepoint.pl tau g jmax : three-leaf values from the ABDKS x-y swap formula (n = 3) after integration by parts,
# compared with Route E's Hodge-side tables engine_h<g>_r3.txt.
#   [x^mu] F_{g,3} = (1/(m1 m2 m3)) [hb^{2g+1}] [z1^{-1}][z2^{-1}][z3^{-1}] ( Q_m1(z1) Q_m2(z2) Q_m3(z3) C ),
#   C = 1/((z1-z2-k12 hb)(z2-z3-k23 hb)(z1-z3+k13 hb)) - 1/((z1-z2+k12 hb)(z2-z3+k23 hb)(z1-z3-k13 hb)),  kij = (mi+mj)/2,
# Laurent expansions at infinity in the region |z1| >> |z2| >> |z3| (z1 first).  Dictionary P2: Xi = sigma^{g-1} [x^mu]F / prod Pi.
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/lib";
use Math::BigRat; use MP;
my ($tau, $g, $jmax, $datafile) = @ARGV; $tau //= 1; $g //= 2; $jmax //= 4;
$datafile //= "$FindBin::Bin/../../E_formulas/engine_h${g}_r3.txt";
$tau = Q($tau); my $e2 = $tau + 1; my $sig = $tau * $e2;
my $H = 2*$g + 1;
sub fact { my $n = shift; my $r = Q(1); $r *= $_ for 2..$n; return $r; }
sub Pi { my $j = shift; my $p = Q(1); $p *= ($j * $tau + $_) for 1..$j-1; return $p / fact($j-1); }
# q_mu(t) coefficients: array over a of arrays over hb-degree
sub qser { my ($mu, $T) = @_; my $n1 = $tau * $mu; my $n2 = $e2 * $mu;
  my @P = map { [ map { Q(0) } 0..$H ] } 0..$T; $P[0][0] = Q(1);
  my $mulf = sub { my ($F) = @_; my @C = map { [ map { Q(0) } 0..$H ] } 0..$T;
    for my $i (0..$T) { for my $k (0..$T-$i) { for my $p (0..$H) { next if $P[$i][$p]->is_zero; for my $q (0..$H-$p) { next if $F->[$k][$q]->is_zero; $C[$i+$k][$p+$q] += $P[$i][$p] * $F->[$k][$q]; } } } }
    @P = @C; };
  for my $kk (1..$n2) { my $m = Q($kk) - ($n2 + 1)/2; my @F = map { [ map { Q(0) } 0..$H ] } 0..$T; $F[0][0] = Q(1); if ($T >= 1) { $F[1][0] = Q(1); $F[1][1] = $m / $e2; } $mulf->(\@F); }
  for my $kk (1..$n1) { my $m = Q($kk) - ($n1 + 1)/2; my @F = map { [ map { Q(0) } 0..$H ] } 0..$T; my $cf = Q(1);
    for my $d (0..$T) { $F[$d][$d] = $cf if $d <= $H; $cf *= -$m / $tau; } $mulf->(\@F); }
  return \@P; }
# variables (z2, z3, hb) -> indices 0,1,2
my $NV = 3;
sub hpoly { my ($arr) = @_; my %p; for my $d (0..$#$arr) { next if $arr->[$d]->is_zero; $p{"0,0,$d"} = $arr->[$d]->copy; } return \%p; }
sub lin { my ($cz2, $cz3, $ch) = @_; my $p = {}; $p = mp_add($p, mp_scal($cz2, mp_var(0,$NV))) if $cz2 != 0; $p = mp_add($p, mp_scal($cz3, mp_var(1,$NV))) if $cz3 != 0; $p = mp_add($p, mp_scal($ch, mp_var(2,$NV))) if $ch != 0; return $p; }
sub trunc_h { my ($p) = @_; my %q; for my $k (keys %$p) { my @e = split /,/, $k; next if $e[2] > $H; $q{$k} = $p->{$k}; } return \%q; }
# complete homogeneous h_p(alpha, beta)
sub hcomp { my ($al, $be, $pmax) = @_; my @apow = (mp_const(1,$NV)); my @bpow = (mp_const(1,$NV));
  for my $i (1..$pmax) { push @apow, trunc_h(mp_mul($apow[-1], $al)); push @bpow, trunc_h(mp_mul($bpow[-1], $be)); }
  my @h; for my $p (0..$pmax) { my $s = {}; for my $a (0..$p) { $s = mp_add($s, trunc_h(mp_mul($apow[$a], $bpow[$p-$a]))); } push @h, $s; } return \@h; }
my %C; open my $fh, '<', $datafile or die "no $datafile";
while (my $line = <$fh>) { chomp $line; my @f = split /\s*\|\s*/, $line; next unless @f >= 4; my ($h, $r, $k, $rest) = @f; next unless $h == $g && $r == 3;
  for my $term (split /\s*;\s*/, $rest) { next unless $term =~ /(\d+),(\d+),(\d+):(\S+)/; $C{$k}{"$1,$2,$3"} = Q($4); } }
close $fh;
sub xi_engine { my ($a1, $a2, $a3) = @_; my $v = Q(0); my $vv = -Q(1) / $sig;
  for my $k (keys %C) { my $s = Q(0); for my $key (keys %{$C{$k}}) { my ($a, $b, $c) = split /,/, $key; $s += $C{$k}{$key} * Q($a1)**$a * Q($a2)**$b * Q($a3)**$c; } $v += $s * $vv**$k; } return $v; }
my ($nb, $nt) = (0, 0);
for my $m1 (1..$jmax) { for my $m2 (1..$m1) { for my $m3 (1..$m2) {
  my @mu = ($m1, $m2, $m3); my $T = 2*($m1 + $m2 + $m3) + 3;
  my @q = map { qser($_, $T) } @mu;
  my ($k12, $k23, $k13) = (Q($m1+$m2)/2, Q($m2+$m3)/2, Q($m1+$m3)/2);
  my $total = Q(0);
  # two terms: sign s = +1 (term B: -k12, -k23, +k13) and s = -1 (term A: +k12, +k23, -k13), weights +1 and -1
  for my $term ([ -1, -1, +1, Q(1) ], [ +1, +1, -1, Q(-1) ]) {
    my ($s12, $s23, $s13, $w) = @$term;
    # z1 step: [z1^{-1}] q1(z1) z1^{m1} / ((z1 - alpha)(z1 - beta)), alpha = z2 + s12 k12 hb, beta = z3 - s13 k13 hb
    my $al = lin(1, 0, -$s12 * $k12); my $be = lin(0, 1, -$s13 * $k13);
    my $hc = hcomp($al, $be, $m1 - 1);
    my $R1 = {};
    for my $p (0..$m1-1) { my $a = $m1 - 1 - $p; $R1 = mp_add($R1, trunc_h(mp_mul(hpoly($q[0][$a]), $hc->[$p]))); }
    # z2 step: [z2^{-1}] R1(z2,z3,hb) q2 z2^{m2} / (z2 - gamma), gamma = z3 + s23 k23 hb  ->  1/(z2-gamma) = Sum_s gamma^s z2^{-s-1}
    my $ga = lin(0, 1, -$s23 * $k23);
    my $emax = 0; for my $key (keys %$R1) { my @e = split /,/, $key; $emax = $e[0] if $e[0] > $emax; }
    my @gpow = (mp_const(1,$NV)); for my $i (1..($emax + $m2 + 1)) { push @gpow, trunc_h(mp_mul($gpow[-1], $ga)); }
    my $R2 = {};
    for my $key (keys %$R1) { my ($e, $f, $hd) = split /,/, $key; my $mono = { "0,$f,$hd" => $R1->{$key} };
      for my $b (0..($e + $m2)) { my $s = $e + $m2 - $b; next if $b > $T; my $qb = hpoly($q[1][$b]); next unless %$qb;
        $R2 = mp_add($R2, trunc_h(mp_mul(mp_mul($mono, $qb), $gpow[$s]))); } }
    # z3 step: [z3^{-1}] R2(z3,hb) q3 z3^{m3}: monomial z3^f needs q3 coefficient c = f + m3 + 1
    for my $key (keys %$R2) { my ($e, $f, $hd) = split /,/, $key; my $c = $f + $m3 + 1; next if $c > $T;
      my $need = $H - $hd; next if $need < 0; my $coef = $q[2][$c][$need]; next if $coef->is_zero;
      $total += $w * $R2->{$key} * $coef; } }
  my $F = $total / ($m1 * $m2 * $m3);
  my $xi_tr = $sig**($g-1) * $F / (Pi($m1) * Pi($m2) * Pi($m3));
  my $xi_h = xi_engine($m1, $m2, $m3);
  $nt++; my $ok = ($xi_tr == $xi_h); $nb++ unless $ok;
  printf "tau=%s g=%d (%d,%d,%d): TR-xy %s  engine %s  %s%s\n", $tau, $g, $m1, $m2, $m3, $xi_tr, $xi_h, $ok ? "OK" : "MISMATCH", ($m1+$m2+$m3 <= 2*$g ? "  [range]" : "");
} } }
print "total $nt, mismatches $nb\n";
