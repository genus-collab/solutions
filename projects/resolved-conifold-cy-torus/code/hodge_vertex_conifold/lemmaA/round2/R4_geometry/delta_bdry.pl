# delta_bdry.pl gmax : the defect Delta_g(j) = B_g(j) - Sum_{g1=0}^{g-1} X_{g-g1}(0) Y_{g1}(j,0)  (Y_g(j,0) = (j-2g+1) X_g(j), Y_0 = 1/j)
#   = (2g-1 - 1/j) X_g(j) - 12 Xi^(1)_g(j) + (1/2)(j-2g+2)(j-2g+3) X_{g-1}(j)      (Prop. M1),
# evaluated with the closed forms X_g = (j-1)_(2g-2)[t^2g] S^2 H^j (Theorem B) and
# Xi^(1)_g = -(1/4)(j-1)_(2g-2)[t^2g] S^2 H^j Q, Q = (1+j/3)(S^2-1) - (j/3) t S S' + t S'/S  (the k=1 part of (6.1)).
# Reports the integer j in 1..3g+1 where Delta_g vanishes; the induction of report s.4 needs Delta_g(j) = 0 for 1 <= j <= 2g-1.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use PJ;
sub Q { Math::BigRat->new(shift) }
my $gmax = shift // 8; $PJ::N = 2*$gmax + 4;
my $f = PJ::logS(); my $tf1 = PJ::stder($f); my $logH = PJ::sadd($f, $tf1, Q(-1));
my $GX = PJ::sexp(PJ::sadd(PJ::sscalp($f, PJ::pc(2)), PJ::sscalp($logH, PJ::pj())));
my $S2 = PJ::sexp(PJ::sscalp($f, PJ::pc(2)));
my $one = PJ::sone();
my $S2m1 = PJ::sadd($S2, $one, Q(-1));
my $tSSp = PJ::sscalp(PJ::stder($S2), PJ::pc(Q(1)/2));        # t S S' = (t/2)(S^2)'
my $Qs = PJ::sadd(PJ::sadd(PJ::sscalp($S2m1, [Q(1), Q(1)/3]), PJ::sscalp($tSSp, [Q(0), Q(-1)/3])), $tf1);
my $GQ = PJ::smul($GX, $Qs);
my %X; for my $g (1..$gmax) { $X{$g} = PJ::pmul(PJ::ffj(-1, 2*$g - 2), $GX->[2*$g]); }
for my $g (2..$gmax) {
  my $X1 = PJ::pscale(PJ::pmul(PJ::ffj(-1, 2*$g - 2), $GQ->[2*$g]), Q(-1)/4);
  my @zeros; my @nz;
  for my $j (1..3*$g + 1) { my $jq = Q($j);
    my $d = (Q(2*$g - 1) - Q(1)/$jq) * PJ::peval($X{$g}, $jq) - 12 * PJ::peval($X1, $jq) + Q(1)/2 * ($jq - 2*$g + 2) * ($jq - 2*$g + 3) * PJ::peval($X{$g-1}, $jq);
    if ($d == 0) { push @zeros, $j } else { push @nz, $j } }
  printf "genus %d: Delta_g(j) = 0 at j = %s;  non-zero at j = %s\n", $g, join(',', @zeros), join(',', @nz); }
