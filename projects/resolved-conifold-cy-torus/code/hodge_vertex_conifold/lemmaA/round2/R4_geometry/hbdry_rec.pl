# hbdry_rec.pl gmax : the boundary-sum identity H(g) of Theorem K1,  B_g(j) = Sum_{g1=0}^{g-1} X_{g-g1}(0) T_{g1}(j)/j,
# tested with B_g computed from Route B's two-leaf REC (recursion (A), as in lam1_rec.pl) -- no closed form for Xi^(1) is used.
# Reports the integers j in 1..3g where it holds (proved for 1 <= j <= 2g-1; j = 2g is the observed extension).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use PJ;
sub Q { Math::BigRat->new(shift) }
$| = 1; my $gmax = shift // 10; $PJ::N = 2*$gmax + 4;
my $f = PJ::logS(); my $tf1 = PJ::stder($f); my $logH = PJ::sadd($f, $tf1, Q(-1));
my $GX = PJ::sexp(PJ::sadd(PJ::sscalp($f, PJ::pc(2)), PJ::sscalp($logH, PJ::pj())));
my $HJ = PJ::sexp(PJ::sscalp($logH, PJ::pj()));                      # H^j
my %X; for my $g (1..$gmax) { $X{$g} = PJ::pmul(PJ::ffj(-1, 2*$g - 2), $GX->[2*$g]); }
my %T; for my $n (0..$gmax) { $T{$n} = PJ::pmul(PJ::ffj(-1, 2*$n), $HJ->[2*$n]); }   # T_n(j) = (j-1)_(2n) [t^2n] H^j
sub c { my $k = shift; return Q(2) / PJ::fact(2*$k) }
sub psi { my ($h1, $k1) = @_; return $h1 == 0 ? PJ::ffj(-1, 2*$k1 - 2) : PJ::pmul(PJ::ffj(-2*$h1 + 1, 2*$k1), $X{$h1}) }
my %B;
for my $h (1..$gmax) {
  my $rhs = PJ::pscale(PJ::pmul([Q(-2*$h + 1), Q(1)], $X{$h}), Q(-(2*$h + 1)));
  for my $k (2..$h) { my $hp = $h + 1 - $k;
    $rhs = PJ::padd($rhs, PJ::pscale(PJ::pmul(PJ::dffj(-2*$hp, 2*$k), PJ::pmul([Q(-2*$hp + 1), Q(1)], $X{$hp})), c($k)));
    $rhs = PJ::padd($rhs, PJ::pscale(PJ::pmul(PJ::ffj(-2*$hp, 2*$k), $B{$hp}), c($k))); }
  $rhs = PJ::padd($rhs, PJ::pscale(PJ::dffj(-1, 2*$h + 1), c($h + 1)));
  for my $h1 (0..$h+1) { for my $k1 (1..$h+1) { for my $h2 (0..$h+1) { my $k2 = $h + 2 - $h1 - $k1 - $h2; next if $k2 < 1;
        my $p = psi($h2, $k2); my $d = $p->[1] // Q(0); next if $d == 0; $rhs = PJ::padd($rhs, PJ::pscale(psi($h1, $k1), c($k1) * c($k2) * $d)); } } }
  my ($q, $rem) = PJ::pdiv_lin($rhs, 2*$h); die "not divisible" if $rem != 0; $B{$h} = PJ::pscale($q, Q(1) / (2*$h + 1));
  my @ok;
  for my $j (1..3*$h) { my $jq = Q($j); my $bs = Q(0); for my $g1 (0..$h-1) { $bs += PJ::peval($X{$h - $g1}, Q(0)) * PJ::peval($T{$g1}, $jq) / $jq; }
    push @ok, $j if PJ::peval($B{$h}, $jq) == $bs; }
  printf "genus %2d: B_g(j) = boundary sum at j = %s\n", $h, join(',', @ok); }
