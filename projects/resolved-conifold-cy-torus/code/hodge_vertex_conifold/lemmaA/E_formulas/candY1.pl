# candY1.pl h r : test candidate Y1 for the pull-back functional
#   sigma^h Y_h(j_1..j_r) =? N! [lambda^{2h}] S1 S2 B^{N+1} Res_w (1+w)^{r-1} prod_i ((1+w)^{a_i} - 1)/a_i (w - (B-1))^{-N-1} dw,
#   a_i = j_i - 1, N = 2h-3+r, B = exp(beta)   (r = 1: this is (6.1)).
# Tested at all integer points 2 <= j_1 <= ... <= j_r <= jmax, as polynomials in sigma.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use YFun; use SerS;
sub R { Math::BigRat->new(shift) }
my ($h, $r, $jmax) = @ARGV; $jmax //= 3*$h;
Hodge::setup($h); my $N = 2*$h - 3 + $r;
my $S12 = SerS::S1S2(); my $bet = SerS::beta(); my $B = SerS::sexp($bet);
my $Bm1 = SerS::sadd($B, SerS::sscale(SerS::sone(), R(-1)));
my $base = SerS::smul($S12, SerS::spow($B, $N + 1));
my @ck; { my $t = $base; for my $k (0..$h) { push @ck, SerS::coeff($t, $h); $t = SerS::smul($t, $Bm1); } }  # [L^h] S1S2 B^{N+1} (B-1)^k
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $x = R(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } return $x; }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
# polynomial in w: coefficients array
sub wmul { my ($p, $q) = @_; my @c = map { R(0) } 0..($#$p + $#$q); for my $i (0..$#$p) { for my $j (0..$#$q) { $c[$i+$j] += $p->[$i] * $q->[$j] } } return \@c; }
sub onepw { my $n = shift; return [ map { binom($n, $_) } 0..$n ]; }
sub cand { my @j = @_; my $poly = onepw($r - 1);
  for my $ji (@j) { my $a = $ji - 1; my $f = onepw($a); $f->[0] -= 1; $f = [ map { $_ / $a } @$f ]; $poly = wmul($poly, $f); }
  my $tot = [R(0)];
  for my $k (0..$h) { my $m = $N + $k; my $c = $m <= $#$poly ? $poly->[$m] : R(0); next if $c == 0;
    $tot = SerS::padd($tot, SerS::pscale($ck[$k], fact($N) * binom($N + $k, $k) * $c)); }
  return $tot; }
sub peq { my ($a, $b) = @_; my $n = $#$a > $#$b ? $#$a : $#$b; for (0..$n) { return 0 if ($a->[$_] // 0) != ($b->[$_] // 0) } return 1; }
sub pstr { my $p = shift; my @t; for my $i (0..$#$p) { next if $p->[$i] == 0; push @t, "($p->[$i])s^$i"; } return @t ? join(' + ', @t) : '0'; }
my @pts = ([]); for (1..$r) { my @n; for my $t (@pts) { my $lo = @$t ? $t->[-1] : 2; for my $x ($lo..$jmax) { push @n, [@$t, $x] } } @pts = @n; }
my ($n, $bad) = (0, 0);
for my $p (@pts) { $n++; my $l = YFun::Ysig(@$p); my $c = cand(@$p);
  unless (peq($l, $c)) { $bad++; print "  R = (@$p): Y = ", pstr($l), "\n           cand = ", pstr($c), "\n" if $bad <= 4; } }
print "genus $h, r = $r: $n points, $bad mismatches\n";
