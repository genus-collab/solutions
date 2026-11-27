package RefXY3;
# Referee 2: Theorem 1 of the R2 report at n = 3 (integer tau >= 1), implemented from scratch by divided differences.
#   [x^mu] F_{g,3} = (1/(m1 m2 m3)) [hb^{2g+1}] [v1^-1][v2^-1][v3^-1] prod_i Q_{m_i}(v_i) * Sum_{3-cycles s} prod_i 1/(v_i - v_s(i) + k_{i s(i)} hb),
#   k_ij = (m_i + m_j)/2, expansions at infinity in |v1| >> |v2| >> |v3|.
# Facts used:  [v^-1] Q(v)/(v - A) = Q_+(A)  (polynomial part at infinity, evaluated at A);
#              [v^-1] Q(v)/((v - A)(v - B)) = (Q_+(A) - Q_+(B))/(A - B) = Sum_a q_a h_{m-1-a}(A, B).
# Polynomials in (v2, v3, hb) are hashes "e2,e3,eh" => BigRat, truncated at hb^{HB}.
use strict; use warnings; use Math::BigRat; use RefR2; use RefXY;
sub Q { Math::BigRat->new(shift) }
my $HB;
sub tadd { my ($p, $q, $c) = @_; $c //= Q(1); for my $k (keys %$q) { my $x = $q->{$k} * $c; if (exists $p->{$k}) { $p->{$k} += $x } else { $p->{$k} = $x } } return $p; }
sub tmul { my ($p, $q) = @_; my %r; for my $a (keys %$p) { my @x = split /,/, $a; for my $b (keys %$q) { my @y = split /,/, $b; my $eh = $x[2] + $y[2]; next if $eh > $HB;
      my $k = ($x[0] + $y[0]) . "," . ($x[1] + $y[1]) . ",$eh"; my $c = $p->{$a} * $q->{$b}; if (exists $r{$k}) { $r{$k} += $c } else { $r{$k} = $c } } } return \%r; }
sub tclean { my $p = shift; for (keys %$p) { delete $p->{$_} if $p->{$_}->is_zero } return $p; }
# linear form c2 v2 + c3 v3 + ch hb
sub lin { my ($c2, $c3, $ch) = @_; my %p; $p{"1,0,0"} = Q($c2) if $c2 != 0; $p{"0,1,0"} = Q($c3) if $c3 != 0; $p{"0,0,1"} = Q($ch) if $ch != 0; return \%p; }
sub powers { my ($L, $n) = @_; my @P = ({ "0,0,0" => Q(1) }); for my $i (1..$n) { push @P, tclean(tmul($P[-1], $L)); } return @P; }
# hb-polynomial (array [p] => coef) as a trivariate polynomial
sub hpoly { my $arr = shift; my %p; for my $d (0..$#$arr) { next if $d > $HB || $arr->[$d]->is_zero; $p{"0,0,$d"} = $arr->[$d]->copy; } return \%p; }
sub F3 { my ($tau, $g, @m) = @_; $HB = 2*$g + 1; my $sumM = $m[0] + $m[1] + $m[2]; my $AM = $sumM + 3;
  my @q = map { RefXY::qcoef($tau, $_, $AM, $HB) } @m;
  my $k = sub { my ($i, $j) = @_; return Q($m[$i] + $m[$j]) / 2; };
  my $total = Q(0);
  # two 3-cycles: s = (1->2->3->1) and s' = (1->3->2->1)  (0-based indices)
  for my $cyc ([1, 2, 0], [2, 0, 1]) {   # cyc->[i] = s(i)
    # factors: for each i, 1/(v_i - v_s(i) + k hb).  v1 (index 0) appears in the factor with i = 0 and in the one with s(i) = 0.
    my $i0 = 0; my $s0 = $cyc->[0]; my ($ip) = grep { $cyc->[$_] == 0 } 0..2;   # ip -> 0
    # factor i=0: 1/(v1 - v_{s0} + k(0,s0) hb) = 1/(v1 - A), A = v_{s0} - k hb
    # factor ip : 1/(v_ip - v1 + k(ip,0) hb) = -1/(v1 - B), B = v_ip + k hb
    my $mkv = sub { my ($idx, $c, $kk) = @_; return lin(($idx == 1 ? 1 : 0), ($idx == 2 ? 1 : 0), $c * $kk); };   # v_idx + c kk hb
    my $A = $mkv->($s0, -1, $k->(0, $s0)); my $B = $mkv->($ip, 1, $k->($ip, 0)); my $sign = Q(-1);
    # divided difference: Sum_a q1_a h_{m1-1-a}(A, B)
    my $m1 = $m[0]; my @PA = powers($A, $m1); my @PB = powers($B, $m1);
    my $DD = {};
    for my $a (0..$m1-1) { my $d = $m1 - 1 - $a; my $h = {}; for my $i (0..$d) { tadd($h, tmul($PA[$i], $PB[$d-$i])); } tclean($h);
      $DD = tadd($DD, tmul(hpoly($q[0][$a]), $h)); }
    tclean($DD);
    # remaining factor (the one with i,s(i) in {1,2}): i = 1 -> s(1) or i = 2 -> s(2)
    my ($ir) = grep { $_ != 0 && $cyc->[$_] != 0 } 0..2; my $sr = $cyc->[$ir];
    # 1/(v_ir - v_sr + k hb); orient with v2 (index 1) as the large variable
    my ($C, $sg2);
    if ($ir == 1) { $C = lin(0, 1, -$k->($ir, $sr)); $sg2 = Q(1); }    # 1/(v2 - (v3 - k hb))
    else { $C = lin(0, 1, $k->($ir, $sr)); $sg2 = Q(-1); }             # 1/(v3 - v2 + k hb) = -1/(v2 - (v3 + k hb))
    # [v2^-1] Q2(v2) DD(v2,...) / (v2 - C) = (Q2 * DD)_+ evaluated at v2 = C
    my %byv2; for my $key (keys %$DD) { my ($e2, $e3, $eh) = split /,/, $key; $byv2{$e2}{"0,$e3,$eh"} = $DD->{$key}; }
    my $m2 = $m[1]; my $maxe2 = 0; $maxe2 = ($_ > $maxe2 ? $_ : $maxe2) for keys %byv2;
    my @PC = powers($C, $m2 + $maxe2);
    my $R2 = {};
    for my $e2 (keys %byv2) { for my $b (0..($m2 + $e2)) { my $pw = $m2 - $b + $e2; next if $pw < 0;   # monomial q2_b v2^{m2-b} * v2^{e2}
        my $t = tmul(hpoly($q[1][$b]), $byv2{$e2}); next unless %$t; $R2 = tadd($R2, tmul($t, $PC[$pw])); } }
    tclean($R2);
    # [v3^-1] Q3(v3) R2(v3, hb)
    my $m3 = $m[2]; my $acc = Q(0);
    for my $key (keys %$R2) { my ($e2, $e3, $eh) = split /,/, $key; die "v2 left" if $e2; my $c = $m3 + $e3 + 1;   # need q3_c with m3 - c + e3 = -1
      next if $c > $AM; my $need = $HB - $eh; next if $need < 0; my $x = $q[2][$c][$need]; next if $x->is_zero; $acc += $R2->{$key} * $x; }
    $total += $sign * $sg2 * $acc; }
  return $total / ($m[0] * $m[1] * $m[2]); }
sub xi3 { my ($tau, $g, @m) = @_; my $sig = $tau * ($tau + 1); my $pp = Q(1); $pp *= RefR2::Pi($_, $tau) for @m; return $sig**($g-1) * F3($tau, $g, @m) / $pp; }
1;
