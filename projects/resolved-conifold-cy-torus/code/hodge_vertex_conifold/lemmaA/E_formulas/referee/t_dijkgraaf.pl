# t_dijkgraaf.pl gmax : Dijkgraaf's two-point function (cited in report 0(d)/(f)(2)) against the referee's WK numbers:
#   Sum_{g} Sum_{a+b=3g-1} <tau_a tau_b>_g x^a y^b = (1/(x+y)) exp((x^3+y^3)/24) Sum_{n>=0} n!/(2n+1)! (x y (x+y)/2)^n   (degree 3g-1 >= 2).
# The right side is computed as a polynomial in each degree after exact division by (x+y).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefWK;
sub Q { Math::BigRat->new(shift) }
my $gmax = shift // 6; my $DEG = 3*$gmax;
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
# bivariate polynomials: hash "a,b" => coef, truncated at total degree DEG+1
sub pmul { my ($p, $q) = @_; my %r; for my $k1 (keys %$p) { my ($a1, $b1) = split /,/, $k1; for my $k2 (keys %$q) { my ($a2, $b2) = split /,/, $k2; next if $a1+$a2+$b1+$b2 > $DEG + 1;
      my $k = ($a1+$a2) . "," . ($b1+$b2); $r{$k} = ($r{$k} // Q(0)) + $p->{$k1} * $q->{$k2}; } } return \%r; }
sub padd { my ($p, $q, $c) = @_; my %r = %$p; for (keys %$q) { $r{$_} = ($r{$_} // Q(0)) + $q->{$_} * $c; } return \%r; }
my $u = { "3,0" => Q('1/24'), "0,3" => Q('1/24') };   # (x^3+y^3)/24
my $E = { "0,0" => Q(1) }; my $t = { "0,0" => Q(1) }; for my $k (1..int(($DEG+1)/3)) { $t = pmul($t, $u); my %s = map { $_ => $t->{$_} / $k } keys %$t; $t = \%s; $E = padd($E, $t, Q(1)); }
my $w = { "2,1" => Q('1/2'), "1,2" => Q('1/2') };     # x y (x+y)/2
my $S = { "0,0" => Q(1) }; my $wp = { "0,0" => Q(1) }; for my $n (1..int(($DEG+1)/3)) { $wp = pmul($wp, $w); $S = padd($S, $wp, fact($n) / fact(2*$n + 1)); }
my $N = pmul($E, $S);   # numerator; divide by (x+y): in each total degree d, c_{a,b} of N = q_{a-1,b} + q_{a,b-1}
my ($n, $bad) = (0, 0);
for my $g (1..$gmax) { my $d = 3*$g - 1;   # quotient degree d, numerator degree d+1
  my @q; for my $a (0..$d) { my $b = $d - $a; my $c = $N->{($a) . "," . ($b+1)} // Q(0); # coefficient of x^a y^{b+1} = q_{a-1,b+1} + q_{a,b}
    $q[$a] = $c - ($a >= 1 ? $q[$a-1] : Q(0)); }
  # check the last relation x^{d+1}: q_{d,0} must equal N_{d+1,0}
  my $last = $N->{($d+1) . ",0"} // Q(0); my $ok = ($q[$d] == $last);
  for my $a (0..$d) { my $v = RefWK::wk($g, $a, $d - $a); $n++; $bad++ if $v != $q[$a]; }
  $bad++ unless $ok; }
printf "Dijkgraaf two-point function vs referee WK numbers: %d coefficients (g = 1..%d), %d mismatches\n", $n, $gmax, $bad;
