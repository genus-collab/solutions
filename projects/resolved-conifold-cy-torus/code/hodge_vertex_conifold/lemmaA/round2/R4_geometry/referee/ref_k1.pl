# ref_k1.pl Gmax : referee's independent re-implementation of the K1 chain (no code from the author's scripts).
#  (1) X_g(j) = (j-1)_(2g-2) [t^2g] S^2 H^j   (Route B, Theorem B),  T_n(j) = (j-1)_(2n) [t^2n] H^j,
#      computed as Sum_k j^k/k! [t^n](A L^k), L = log H  (exponential-in-j expansion, not a series of j-polynomials).
#  (2) the FULL bivariate two-leaf Hurwitz polynomial Y_h(j,y) from Route B's REC for R = (j,y) (genus-0 terms in polynomial form),
#      by exact division by (j + y - 2h); the two-block part is computed as the explicit double sum AND compared with
#      Sum_{m1+m2=h} T_{m1}(j) T_{m2}(y) (step 4 of Prop. M2).
#  (3) consistency: Y_h(j,0) = (j-2h+1) X_h(j) = T_h(j)/j (dilaton), Y_h(j,1) = j X_h(j) (string)  -- these are NOT used as input.
#  (4) B_h = [y^1] Y_h; the author's recursion (A) re-checked as a polynomial identity;  M3: B_h(j) = BS_h(j) at integers.
#  (5) Xi^(1)_h(j) from Mumford (Prop. M1); K1 in the range; K1' at j = 2h-1; M3' at j = 2h; the closed form C1.
#  (6) comparison with tables given on the command line through env TABLES (dir:pattern), done in a separate script.
# Output: referee/out/ref_k1_<Gmax>.txt (summary) and referee/out/xi1_ref.txt, Y_ref_h<h>.txt (polynomials).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RL; use Time::HiRes qw(time);
sub Q { Math::BigRat->new(shift) }
$| = 1; my $G = shift // 8; my $t0 = time;
$RL::NT = 2*$G + 6;
my $L = RL::logH(); my $f = RL::logS(); my $one = RL::szero(); $one->[0] = Q(1);
my $S2 = RL::sexp(RL::ssc($f, Q(2)));
my (%X, %T);
for my $g (1..$G+1) { $X{$g} = RL::pmul(RL::ffj(-1, 2*$g - 2), RL::coeff_expj($S2, $L, 2*$g)) }
for my $n (0..$G+1) { $T{$n} = RL::pmul(RL::ffj(-1, 2*$n), RL::coeff_expj($one, $L, 2*$n)) }
my %chi; $chi{$_} = -RL::bern(2*$_) / (2*$_) for 1..$G+1;
my %zeta; $zeta{$_} = (2*$_ - 1) * RL::bern(2*$_) / (2*$_) for 1..$G+1;
my @msg;
# sanity: X_g(0) = chi_g (Harer-Zagier), T_n(0) = 0, T_n'(0) = zeta_n
for my $g (1..$G) { push @msg, "X_$g(0) != chi_$g" if RL::pev($X{$g}, Q(0)) != $chi{$g}; push @msg, "T_$g(0) != 0" if RL::pev($T{$g}, Q(0)) != 0; push @msg, "T_$g'(0) != zeta_$g" if (RL::pder($T{$g})->[0]) != $zeta{$g}; }
sub c { my $k = shift; return Q(2) / RL::fact(2*$k) }
# psi_{h,k}(x) = (x-2h+1)_(2k) X_h(x) for h >= 1;  (x-1)_(2k-2) for h = 0
sub psiu { my ($h, $k) = @_; return $h == 0 ? RL::ffj(-1, 2*$k - 2) : RL::pmul(RL::ffj(-2*$h + 1, 2*$k), $X{$h}) }
my %Y;   # bivariate, h >= 1
my ($bad_tb, $bad_cons) = (0, 0);
for my $h (1..$G) {
  my $rhs = [[Q(0)]];
  for my $k (2..$h) { my $hp = $h + 1 - $k; $rhs = RL::badd($rhs, RL::bsc(RL::bmul(RL::ffjy(-2*$hp, 2*$k), $Y{$hp}), c($k))) }
  $rhs = RL::badd($rhs, RL::bsc(RL::ffjy(-1, 2*$h + 1), c($h + 1)));            # genus 0: (j+y)_(2h+2)/(j+y)
  my $two = [[Q(0)]];
  for my $h1 (0..$h+1) { for my $k1 (1..$h+1) { for my $h2 (0..$h+1) { my $k2 = $h + 2 - $h1 - $k1 - $h2; next if $k2 < 1;
        my $pj = psiu($h1, $k1); my $py = psiu($h2, $k2); my $b = [ map { RL::psc($pj, $_) } @$py ];   # P(j) Q(y)
        $two = RL::badd($two, RL::bsc(RL::btrim($b), c($k1) * c($k2))); } } }
  # check: two-block part = Sum_{m1+m2=h} T_{m1}(j) T_{m2}(y)
  my $tt = [[Q(0)]]; for my $m1 (0..$h) { my $m2 = $h - $m1; $tt = RL::badd($tt, RL::btrim([ map { RL::psc($T{$m1}, $_) } @{$T{$m2}} ])) }
  my $dd = RL::badd($two, $tt, Q(-1)); for my $p (@$dd) { $bad_tb++ unless RL::pzero($p) }
  $rhs = RL::badd($rhs, $two);
  my ($q, $rem) = RL::bdivlin($rhs, 2*$h); die "genus $h: REC right side not divisible by (j+y-2h)" unless RL::pzero($rem);
  $Y{$h} = RL::bsc($q, Q(1) / (2*$h + 1));
  # consistency (string/dilaton) -- not used as input
  my $y0 = $Y{$h}[0]; my $y1 = [Q(0)]; for my $b (0..$#{$Y{$h}}) { $y1 = RL::padd($y1, $Y{$h}[$b]) }
  $bad_cons++ unless RL::peq($y0, RL::pmul([Q(-2*$h + 1), Q(1)], $X{$h}));
  $bad_cons++ unless RL::peq(RL::pmul($y0, [Q(0), Q(1)]), $T{$h});
  $bad_cons++ unless RL::peq($y1, RL::pmul([Q(0), Q(1)], $X{$h}));
  printf "  genus %2d: Y_h computed (deg_y %d) [%.1f s]\n", $h, $#{$Y{$h}}, time - $t0; }
push @msg, "two-block simplification failed in $bad_tb y-coefficients" if $bad_tb;
push @msg, "string/dilaton consistency failed $bad_cons times" if $bad_cons;
print "two-block part = Sum T_m1(j) T_m2(y): ", ($bad_tb ? "FAILS" : "holds"), " for h <= $G;  Y_h(j,0) = (j-2h+1)X_h = T_h/j and Y_h(j,1) = j X_h: ", ($bad_cons ? "FAIL" : "hold"), "\n";
# B_h and recursion (A) as a polynomial identity
my %B; $B{$_} = $Y{$_}[1] // [Q(0)] for 1..$G;
my $badA = 0;
for my $h (1..$G) {
  my $lhs = RL::psc(RL::padd(RL::pmul([Q(-2*$h + 1), Q(1)], $X{$h}), RL::pmul([Q(-2*$h), Q(1)], $B{$h})), Q(2*$h + 1));
  my $r = [Q(0)];
  for my $hp (1..$h-1) { my $n = 2*$h + 2 - 2*$hp; my ($tq, $tr) = RL::pdivlin($T{$hp}, Q(0)); die unless $tr == 0;
    $r = RL::padd($r, RL::psc(RL::padd(RL::pmul(RL::pder(RL::ffj(-2*$hp, $n)), $tq), RL::pmul(RL::ffj(-2*$hp, $n), $B{$hp})), c($h + 1 - $hp))); }
  # h'' = 0: c_{h+1} [ (x)_(2h+2)'|_j / j - (j)_(2h+2)/j^2 ] = c_{h+1} d/dj (j-1)_(2h+1)
  $r = RL::padd($r, RL::psc(RL::pder(RL::ffj(-1, 2*$h + 1)), c($h + 1)));
  for my $m (1..$h) { $r = RL::padd($r, RL::psc($T{$h - $m}, $zeta{$m})) }
  $badA++ unless RL::peq($lhs, $r); }
print "recursion (A) of Prop. M2 as a polynomial identity in j, h = 1..$G: ", ($badA ? "$badA FAILURES" : "holds"), "\n";
# M3 and Mumford reduction
my %X1;
for my $h (1..$G) {
  my ($q0, $r0) = RL::pdivlin(RL::padd($X{$h}, [ -$chi{$h} ]), Q(0)); die unless $r0 == 0;          # (X_h(j) - X_h(0))/j
  my $t = RL::padd(RL::padd(RL::psc($X{$h}, Q(2*$h - 1)), $B{$h}, Q(-1)), $q0, Q(-1));
  if ($h >= 2) { $t = RL::padd($t, RL::psc(RL::pmul(RL::pmul([Q(-2*$h + 2), Q(1)], [Q(-2*$h + 3), Q(1)]), $X{$h-1}), Q(1)/2)) } else { $t = RL::padd($t, [Q(1)/2]) }
  for my $g1 (1..$h-1) { $t = RL::padd($t, RL::psc(RL::pmul([Q(-2*$g1 + 1), Q(1)], $X{$g1}), $chi{$h - $g1})) }
  $X1{$h} = RL::psc($t, Q(1)/12); }
# closed form C1: -(1/4)(j-1)_(2h-2) [t^2h] S^2 H^j [(1+j/3)(S^2-1) - (j/3) t S S' + t S'/S]
my $S2m1 = RL::sadd($S2, $one, Q(-1)); my $tSSp = RL::ssc(RL::stdt($S2), Q(1)/2); my $tlogSp = RL::stdt($f);
my $A0 = RL::smul($S2, RL::sadd($S2m1, $tlogSp));                       # j^0 part of S^2 [..]
my $A1 = RL::smul($S2, RL::ssc(RL::sadd($S2m1, $tSSp, Q(-1)), Q(1)/3));    # coefficient of j in [..], times S^2
open my $o, '>', "$FindBin::Bin/out/xi1_ref.txt" or die;
my @sum;
for my $h (1..$G) {
  my $cf = RL::padd(RL::coeff_expj($A0, $L, 2*$h), RL::pmul([Q(0), Q(1)], RL::coeff_expj($A1, $L, 2*$h)));
  $cf = RL::psc(RL::pmul(RL::ffj(-1, 2*$h - 2), $cf), Q(-1)/4);
  my $c1 = RL::peq($cf, $X1{$h});
  my @nzr = grep { RL::pev($X1{$h}, Q($_)) != 0 } 2..(2*$h - 2);
  my $k1p = 1; if ($h >= 2) { my $jj = Q(2*$h - 1); $k1p = (12 * RL::pev($X1{$h}, $jj) == (Q(2*$h - 1) - 1/$jj) * RL::pev($X{$h}, $jj) + RL::pev($X{$h-1}, $jj)) ? 1 : 0 }
  my @m3; my @m3x; for my $j (1..2*$h + 2) { my $jq = Q($j); my $bs = Q(0); for my $g1 (0..$h-1) { $bs += $chi{$h - $g1} * RL::pev($T{$g1}, $jq) / $jq }
    if (RL::pev($B{$h}, $jq) == $bs) { push @m3, $j } else { push @m3x, $j } }
  my $m3ok = (join(',', @m3) =~ /^1(,\d+)*$/ && scalar(grep { $_ <= 2*$h - 1 } @m3) == 2*$h - 1) ? 1 : 0;
  my $line = sprintf "genus %2d: C1 %s; K1 (j=2..%d) %s; K1' %s; M3 holds at j = %s (fails at %s)", $h, ($c1 ? "equal" : "DIFFERS"), 2*$h-2, (@nzr ? "FAILS at @nzr" : "zero"), ($k1p ? "ok" : "FAILS"), join(',', @m3), join(',', @m3x);
  print "$line\n"; push @sum, $line;
  print $o "$h | 1 | 1 | ", RL::pstr($X1{$h}), "\n"; }
close $o;
for my $h (1..$G) { open my $oy, '>', "$FindBin::Bin/out/Y_ref_h${h}.txt" or die; my @terms; for my $b (0..$#{$Y{$h}}) { my $p = $Y{$h}[$b]; for my $a (0..$#$p) { push @terms, "$a,$b:$p->[$a]" if $p->[$a] != 0 } }
  print $oy "$h | 2 | 0 | ", join(' ; ', @terms), "\n"; close $oy;
  open my $ox, '>', "$FindBin::Bin/out/X_ref_h${h}.txt" or die; print $ox "$h | 1 | 0 | ", RL::pstr($X{$h}), "\n"; close $ox; }
print "sanity messages: ", (@msg ? join('; ', @msg) : 'none'), "\n";
printf "done [%.1f s]\n", time - $t0;
