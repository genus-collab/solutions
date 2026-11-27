# ref_rec_A.pl Hmax : B_h from the y-derivative recursion (A) written in the T_n / zeta form of the report
# (the author's scripts use the explicit two-block double sum instead); then Prop. M1, K1, K1', M3, M3' and the closed form C1
# for h <= Hmax.  (A) itself was checked as a polynomial identity against the full bivariate REC in ref_k1.pl.
#   (2h+1)[(j-2h+1) X_h + (j-2h) B_h] = Sum_{h''=1}^{h-1} c_{h+1-h''} [ d/dx (x)_(n)|_{x=j-2h''} T_{h''}/j + (j-2h'')_(n) B_{h''} ]
#        + c_{h+1} d/dj (j-1)_(2h+1) + Sum_{m=1}^h zeta_m T_{h-m},     n = 2h+2-2h''.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RL; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my $H = shift // 16; my $t0 = time;
$RL::NT = 2*$H + 4;
my $L = RL::logH(); my $f = RL::logS(); my $one = RL::szero(); $one->[0] = Q(1); my $S2 = RL::sexp(RL::ssc($f, Q(2)));
my (%X, %T, %chi, %zeta);
for my $g (1..$H) { $X{$g} = RL::pmul(RL::ffj(-1, 2*$g - 2), RL::coeff_expj($S2, $L, 2*$g)) }
for my $n (0..$H) { $T{$n} = RL::pmul(RL::ffj(-1, 2*$n), RL::coeff_expj($one, $L, 2*$n)) }
for my $m (1..$H) { $chi{$m} = -RL::bern(2*$m) / (2*$m); $zeta{$m} = (2*$m - 1) * RL::bern(2*$m) / (2*$m) }
printf "series and X, T ready [%.1f s]\n", time - $t0;
sub c { my $k = shift; return Q(2) / RL::fact(2*$k) }
my %B;
for my $h (1..$H) {
  my $r = [Q(0)];
  for my $hp (1..$h-1) { my $n = 2*$h + 2 - 2*$hp; my ($tq, $tr) = RL::pdivlin($T{$hp}, Q(0)); die unless $tr == 0;
    $r = RL::padd($r, RL::psc(RL::padd(RL::pmul(RL::pder(RL::ffj(-2*$hp, $n)), $tq), RL::pmul(RL::ffj(-2*$hp, $n), $B{$hp})), c($h + 1 - $hp))); }
  $r = RL::padd($r, RL::psc(RL::pder(RL::ffj(-1, 2*$h + 1)), c($h + 1)));
  for my $m (1..$h) { $r = RL::padd($r, RL::psc($T{$h - $m}, $zeta{$m})) }
  $r = RL::padd($r, RL::pmul([Q(-2*$h + 1), Q(1)], $X{$h}), Q(-(2*$h + 1)));
  my ($q, $rem) = RL::pdivlin($r, Q(2*$h)); die "genus $h: not divisible by (j-2h), remainder $rem" if $rem != 0;
  $B{$h} = RL::psc($q, Q(1) / (2*$h + 1)); }
printf "B_h ready [%.1f s]\n", time - $t0;
my $S2m1 = RL::sadd($S2, $one, Q(-1)); my $tSSp = RL::ssc(RL::stdt($S2), Q(1)/2); my $tlogSp = RL::stdt($f);
my $A0 = RL::smul($S2, RL::sadd($S2m1, $tlogSp)); my $A1 = RL::smul($S2, RL::ssc(RL::sadd($S2m1, $tSSp, Q(-1)), Q(1)/3));
my ($allC1, $allK1, $allK1p, $allM3, $allM3p) = (1, 1, 1, 1, 1);
for my $h (1..$H) {
  my ($q0, $r0) = RL::pdivlin(RL::padd($X{$h}, [ -$chi{$h} ]), Q(0)); die unless $r0 == 0;
  my $t = RL::padd(RL::padd(RL::psc($X{$h}, Q(2*$h - 1)), $B{$h}, Q(-1)), $q0, Q(-1));
  if ($h >= 2) { $t = RL::padd($t, RL::psc(RL::pmul(RL::pmul([Q(-2*$h + 2), Q(1)], [Q(-2*$h + 3), Q(1)]), $X{$h-1}), Q(1)/2)) } else { $t = RL::padd($t, [Q(1)/2]) }
  for my $g1 (1..$h-1) { $t = RL::padd($t, RL::psc(RL::pmul([Q(-2*$g1 + 1), Q(1)], $X{$g1}), $chi{$h - $g1})) }
  my $X1 = RL::psc($t, Q(1)/12);
  my $cf = RL::padd(RL::coeff_expj($A0, $L, 2*$h), RL::pmul([Q(0), Q(1)], RL::coeff_expj($A1, $L, 2*$h)));
  $cf = RL::psc(RL::pmul(RL::ffj(-1, 2*$h - 2), $cf), Q(-1)/4);
  my $c1 = RL::peq($cf, $X1); $allC1 &&= $c1;
  my @nz = grep { RL::pev($X1, Q($_)) != 0 } 2..(2*$h - 2); $allK1 &&= !@nz;
  my $k1p = 1; if ($h >= 2) { my $jj = Q(2*$h - 1); $k1p = (12 * RL::pev($X1, $jj) == (Q(2*$h - 1) - 1/$jj) * RL::pev($X{$h}, $jj) + RL::pev($X{$h-1}, $jj)) ? 1 : 0 } $allK1p &&= $k1p;
  my @ok; my @no; for my $j (1..2*$h + 1) { my $jq = Q($j); my $bs = Q(0); for my $g1 (0..$h-1) { $bs += $chi{$h - $g1} * RL::pev($T{$g1}, $jq) / $jq }
    if (RL::pev($B{$h}, $jq) == $bs) { push @ok, $j } else { push @no, $j } }
  my $m3 = (scalar(grep { $_ <= 2*$h - 1 } @ok) == 2*$h - 1) ? 1 : 0; $allM3 &&= $m3; my $m3p = (grep { $_ == 2*$h } @ok) ? 1 : 0; $allM3p &&= $m3p;
  my $m3pp = 1; if ($h >= 2) { my $jj = Q(2*$h); $m3pp = (12 * RL::pev($X1, $jj) == 3 * RL::pev($X{$h-1}, $jj)) ? 1 : 0 }
  printf "genus %2d: C1 %s | K1 %s | K1' %s | M3 on 1..%d %s | M3' (j=2h) %s, 12Xi1(2h)=3X_{h-1}(2h) %s | M3 fails at %s [%.1f s]\n", $h, ($c1 ? 'ok' : 'FAIL'), (@nz ? "FAIL @nz" : 'ok'), ($k1p ? 'ok' : 'FAIL'),
    2*$h - 1, ($m3 ? 'ok' : 'FAIL'), ($m3p ? 'ok' : 'no'), ($m3pp ? 'ok' : 'no'), join(',', @no), time - $t0; }
printf "SUMMARY h <= %d: C1 %s, K1 %s, K1' %s, M3 %s, M3' %s\n", $H, map { $_ ? 'all ok' : 'FAILURES' } ($allC1, $allK1, $allK1p, $allM3, $allM3p);
