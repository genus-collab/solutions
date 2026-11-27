# residual_rank.pl h sigma : r = 2 residual  Delta(j1,j2) = sigma^h Y_h(j1,j2) - Y1(j1,j2)  at a numeric sigma;
# prints the coefficient matrix of Delta in the monomials (j1-1)^a (j2-1)^b and its rank (exact Gaussian elimination).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use YFun; use SerS;
sub R { Math::BigRat->new(shift) }
my ($h, $sg) = @ARGV; $sg = R($sg // 2);
Hodge::setup($h); my $r = 2; my $N = 2*$h - 3 + $r; my $D = 3*$h - 3;
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
sub binom { my ($n, $k) = @_; return R(0) if $k < 0 || $k > $n; my $x = R(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } return $x; }
# Y as a polynomial in (j1, j2): coefficient of j1^(a-1) j2^(b-1) is Ypsi(a,b) evaluated at sigma
my %Y; for my $a (1..$D+1) { for my $b (1..$D+2-$a) { my $v = YFun::Ypsi($a, $b); my $s = R(0);
    for my $k (0..$#$v) { $s += (($k % 2) ? -1 : 1) * $v->[$k] * $sg ** ($h - $k); } $Y{($a-1).",".($b-1)} = $s; } }
# candidate at numeric sigma, as polynomial in a1 = j1-1, a2 = j2-1 : N! Sum_k C(N+k,k) [L^h](S1S2 B^{N+1}(B-1)^k) c_{N+k}(a1,a2)
my $S12 = SerS::S1S2(); my $B = SerS::sexp(SerS::beta()); my $Bm1 = SerS::sadd($B, SerS::sscale(SerS::sone(), R(-1)));
my @ck; { my $t = SerS::smul($S12, SerS::spow($B, $N + 1)); for my $k (0..$h) { push @ck, SerS::peval(SerS::coeff($t, $h), $sg); $t = SerS::smul($t, $Bm1); } }
# c_m(a1,a2) = [w^m] (1+w) f(a1) f(a2), f(a) = Sum_{m>=1} C(a,m)/a w^m ; C(a,m)/a as polynomial in a (coeff list)
sub polymul { my ($p, $q) = @_; my @c = map { R(0) } 0..($#$p + $#$q); for my $i (0..$#$p) { for my $j (0..$#$q) { $c[$i+$j] += $p->[$i] * $q->[$j] } } return \@c; }
sub Cm_over_a { my $m = shift; my $p = [R(1)]; for my $i (1..$m-1) { $p = polymul($p, [R(-$i), R(1)]); } return [ map { $_ / fact($m) } @$p ]; }  # (a-1)...(a-m+1)/m!
my %cand;  # key "e1,e2" (powers of a1, a2)
for my $k (0..$h) { my $m = $N + $k; my $w = fact($N) * binom($N + $k, $k) * $ck[$k];
  for my $m1 (1..$m) { for my $eps (0, 1) { my $m2 = $m - $m1 - $eps; next if $m2 < 1;   # (1+w) = w^0 + w^1
      my $p1 = Cm_over_a($m1); my $p2 = Cm_over_a($m2);
      for my $e1 (0..$#$p1) { for my $e2 (0..$#$p2) { my $c = $p1->[$e1] * $p2->[$e2] * $w; next if $c == 0; $cand{"$e1,$e2"} = ($cand{"$e1,$e2"} // R(0)) + $c; } } } } }
# Y in powers of (j-1): re-expand j^e = Sum C(e,i) (j-1)^i
my %Ya; for my $key (keys %Y) { my ($e1, $e2) = split /,/, $key; for my $i1 (0..$e1) { for my $i2 (0..$e2) { my $c = $Y{$key} * binom($e1, $i1) * binom($e2, $i2); $Ya{"$i1,$i2"} = ($Ya{"$i1,$i2"} // R(0)) + $c; } } }
my %Dl; for my $k (keys %Ya, keys %cand) { $Dl{$k} = ($Ya{$k} // R(0)) - ($cand{$k} // R(0)); }
my $M = 0; for (keys %Dl) { my ($a, $b) = split /,/; $M = $a if $a > $M; $M = $b if $b > $M; }
my @mat; for my $a (0..$M) { push @mat, [ map { $Dl{"$a,$_"} // R(0) } 0..$M ]; }
print "genus $h, sigma = $sg: residual Delta in powers (j1-1)^a (j2-1)^b, a,b = 0..$M:\n";
for my $a (0..$M) { print "  ", join(' ', map { sprintf("%s", $_) } @{$mat[$a]}), "\n"; }
# rank
my @A = map { [ map { $_->copy } @$_ ] } @mat; my $rank = 0; my $nc = $M + 1;
for my $c (0..$nc-1) { my $piv; for my $i ($rank..$#A) { if ($A[$i][$c] != 0) { $piv = $i; last } } next unless defined $piv;
  @A[$rank, $piv] = @A[$piv, $rank]; for my $i (0..$#A) { next if $i == $rank || $A[$i][$c] == 0; my $f = $A[$i][$c] / $A[$rank][$c]; $A[$i][$_] -= $f * $A[$rank][$_] for 0..$nc-1; } $rank++; }
print "rank of the residual matrix: $rank\n";
