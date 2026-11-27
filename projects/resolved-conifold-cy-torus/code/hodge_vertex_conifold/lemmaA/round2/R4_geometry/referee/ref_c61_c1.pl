# ref_c61_c1.pl Hmax : is C1 exactly the order-1/sigma part of (6.1)?  For h <= Hmax, extract from the right side of (6.1),
#   sigma^h Xi_h(j) = (j-1)_(2h-2) [l^{2h}] exp( Sum_m c_m l^{2m} (p_{2m} - j sigma p_{2m-1}) ),
# the coefficients of sigma^h (-> Xi^(0)) and sigma^(h-1) (-> -Xi^(1)) as polynomials in j, and compare with Theorem B and with
#   C1:  Xi^(1)_h(j) = -(1/4)(j-1)_(2h-2) [t^2h] S^2 H^j [(1+j/3)(S^2-1) - (j/3) t S S' + t S'/S].
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RL;
sub Q { Math::BigRat->new(shift) }
my $H = shift // 12; $RL::NT = 2*$H + 4;
my $f = RL::logS(); my $L = RL::logH(); my $one = RL::szero(); $one->[0] = Q(1); my $S2 = RL::sexp(RL::ssc($f, Q(2)));
my $S2m1 = RL::sadd($S2, $one, Q(-1)); my $tSSp = RL::ssc(RL::stdt($S2), Q(1)/2); my $tlogSp = RL::stdt($f);
my $A0 = RL::smul($S2, RL::sadd($S2m1, $tlogSp)); my $A1 = RL::smul($S2, RL::ssc(RL::sadd($S2m1, $tSSp, Q(-1)), Q(1)/3));
my @pn = ([Q(2)], [Q(1)]); for my $n (2..2*$H) { $pn[$n] = RL::padd($pn[$n-1], [Q(0), @{$pn[$n-2]}]) }
my @e; for my $m (1..$H) { my $cm = -RL::bern(2*$m) / (Q(2*$m) * RL::fact(2*$m));
  $e[$m] = RL::badd([ map { [$_ * $cm] } @{$pn[2*$m]} ], [ [Q(0)], map { [Q(0), -$_ * $cm] } @{$pn[2*$m - 1]} ]) }
my @a = ([[Q(1)]]); for my $n (1..$H) { my $s = [[Q(0)]]; for my $m (1..$n) { $s = RL::badd($s, RL::bsc(RL::bmul($e[$m], $a[$n - $m]), Q($m))) } $a[$n] = RL::bsc($s, Q(1) / $n) }
for my $h (1..$H) { my $rhs = RL::bmulj($a[$h], RL::ffj(-1, 2*$h - 2));
  my $x0 = $rhs->[$h] // [Q(0)]; my $x1 = RL::psc($rhs->[$h - 1] // [Q(0)], Q(-1));
  my $Xh = RL::pmul(RL::ffj(-1, 2*$h - 2), RL::coeff_expj($S2, $L, 2*$h));
  my $cf = RL::psc(RL::pmul(RL::ffj(-1, 2*$h - 2), RL::padd(RL::coeff_expj($A0, $L, 2*$h), RL::pmul([Q(0), Q(1)], RL::coeff_expj($A1, $L, 2*$h)))), Q(-1)/4);
  my $deg = $#$rhs;
  printf "genus %2d: (6.1) sigma-degree %d; sigma^h coefficient = Theorem B: %s; sigma^(h-1) coefficient = -C1: %s\n", $h, $deg, (RL::peq($x0, $Xh) ? 'yes' : 'NO'), (RL::peq($x1, $cf) ? 'yes' : 'NO'); }
