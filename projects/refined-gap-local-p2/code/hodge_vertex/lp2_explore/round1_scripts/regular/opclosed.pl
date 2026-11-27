use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RegLib; use Math::BigRat;
# Test 7, part D: closed forms of the WKB recursion operators of the T_F PF equation in s = 1/(L+alpha), for n = 1..NMAX:
#   Op1 = -s^-1 d/ds + s^-2                       (i.e. Op1(s sigma) = -sigma')
#   Op2(s^n) = -(n-1)(2n-1)/(3 alpha) s^(n-2) + (3n^2-3n+1)/9 s^(n-1)
#   Op3(s^n) = -n(n-1)^2/(9 alpha^2) s^(n-2) + n(3n^2-3n+1)/(27 alpha) s^(n-1) - n^3/27 s^n
# Each s-coefficient of s^(3-n) q^3 Op_j(s^n) is a polynomial in n of degree <= j (theta applied j times to s^n), and the
# s-range is bounded, so agreement for NMAX >= 4 values of n proves the closed forms for all n >= 1 (here NMAX = 15).
my ($E1,$E2,$NMAX) = @ARGV; $E1=QQ($E1); $E2=QQ($E2); $NMAX //= 15; my $E=$E1+$E2; my $al=$E/3;
my $a = {-1 => QQ(1)}; my $Y = {-1 => QQ(1), 0 => -$al}; my $q = {0 => QQ(1), 1 => -$al}; my $q3 = lp_mul($q, lp_mul($q,$q));
sub th { lp_theta($_[0],$al) }
sub q3Op1 { my $f=shift; my $fp = th($f);
  my $A = lp_add(lp_scale(lp_mul(lp_mul($a,$a),$fp),QQ(3)), lp_scale(lp_mul(lp_mul($a,th($a)),$f),QQ(3)));
  my $B = lp_add(lp_scale(lp_mul(lp_mul($Y,$Y),$fp),QQ(3)), lp_scale(lp_mul(lp_mul($Y,th($Y)),$f),QQ(3)), lp_mul(lp_mul($Y,$Y),$f));
  return lp_add(lp_mul($q3,$A), lp_scale($B,QQ(-1))); }
sub q3Op2 { my $f=shift; my $fp = th($f); my $fpp = th($fp);
  my $A = lp_add(lp_scale(lp_mul($a,$fpp),QQ(3)), lp_scale(lp_mul(th($a),$fp),QQ(3)), lp_mul(th(th($a)),$f));
  my $B = lp_add(lp_scale(lp_mul($Y,$fpp),QQ(3)), lp_scale(lp_mul(th($Y),$fp),QQ(3)), lp_mul(th(th($Y)),$f),
                 lp_scale(lp_mul($Y,$fp),QQ(2)), lp_mul(th($Y),$f), lp_scale(lp_mul($Y,$f),QQ(2)/9));
  return lp_add(lp_mul($q3,$A), lp_scale($B,QQ(-1))); }
sub q3Op3 { my $f=shift; my $f1 = th($f); my $f2 = th($f1); my $f3 = th($f2);
  return lp_add(lp_mul($q3,$f3), lp_scale(lp_add($f3,$f2,lp_scale($f1,QQ(2)/9)),QQ(-1))); }
my $bad = 0;
for my $n (1..$NMAX) {
  my $N = QQ($n);
  my $o1 = lp_add(q3Op1({$n => QQ(1)}), lp_scale(lp_mul($q3, {$n-2 => -$N, }), QQ(1)), lp_scale(lp_mul($q3,{$n-2 => QQ(1)}),QQ(1)));
  # Op1(s^n) = -n s^(n-2) + s^(n-2) = (1-n) s^(n-2); so q3Op1(s^n) - q3*(1-n)s^(n-2) should vanish: o1 = q3Op1 + q3*(n-1)... (built above as q3Op1 - n q3 s^(n-2) + q3 s^(n-2))
  my $o1c = lp_add(q3Op1({$n => QQ(1)}), lp_scale(lp_mul($q3,{$n-2 => QQ(1)}), $N-1));
  my $c2 = { $n-2 => -($N-1)*(2*$N-1)/(3*$al), $n-1 => (3*$N*$N-3*$N+1)/9 };
  my $c3 = { $n-2 => -$N*($N-1)**2/(9*$al*$al), $n-1 => $N*(3*$N*$N-3*$N+1)/(27*$al), $n => -$N**3/27 };
  my $d2 = lp_add(q3Op2({$n => QQ(1)}), lp_scale(lp_mul($q3, lp_clean($c2)),QQ(-1)));
  my $d3 = lp_add(q3Op3({$n => QQ(1)}), lp_scale(lp_mul($q3, lp_clean($c3)),QQ(-1)));
  my $r = sprintf "n=%2d: Op1 closed form %s; Op2 closed form %s; Op3 closed form %s", $n, ((keys %$o1c)?'FAIL':'ok'), ((keys %$d2)?'FAIL':'ok'), ((keys %$d3)?'FAIL':'ok');
  $bad++ if (keys %$o1c) || (keys %$d2) || (keys %$d3);
  print "$r\n";
}
print "closed forms: ", ($bad ? "$bad FAILURES" : "all exact for n = 1..$NMAX"), "  (e1,e2) = ($E1,$E2)\n";
