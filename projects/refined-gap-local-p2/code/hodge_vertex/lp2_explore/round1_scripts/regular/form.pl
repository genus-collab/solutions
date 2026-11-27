use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RegLib; use Math::BigRat;
# Test 7 (regular colours), part A: the form cPP_k = Sqrt[3]^p_k sq^-1 N_k(L)/(E+3L)^(k+1), deg N_k <= k,
# rewritten in s := 1/(L + alpha) = 3/(E + 3L)   (alpha = E/3).  Form <=> S_k(s) := cPP_k |_{sq=1, Sqrt3 -> 1} in s Q[s], deg <= k+1.
# Output: textual denominator exponents and max L-power, the Sqrt[3] power, exact s-coefficients of S_k, and checks.
# Writes the exact S_k (with the convention cPP_k = Sqrt[3] sq^-1 * Shat_k, Shat_k rational) to Sk_<e1>_<e2>.txt for the other scripts.
my ($f,$E1,$E2) = @ARGV; $E1 = QQ($E1); $E2 = QQ($E2); my $E = $E1+$E2; my $al = $E/3;
my @parts = parse_cpp($f); my $N = $#parts;
print "file $f: ", $N+1, " coefficients z^0..z^$N; (e1,e2) = ($E1,$E2)\n";
my $bs = chr(92);
open my $out, '>', "$FindBin::Bin/Sk_${E1}_${E2}.txt" or die;
for my $k (0..$N) {
  my $t = $parts[$k]; $t =~ s/\Q$bs\E//g;
  my $pat = '([Epsilon][1]+[Epsilon][2]+3*LLL[0][y])';
  my @m; my $tt = $t; while ($tt =~ s/\Q$pat\E(\^(\d+))?/#/) { push @m, (defined $2 ? $2 : 1) }
  my $maxL = 0; while ($tt =~ /LLL\[0\]\[y\](\^(\d+))?/g) { my $d = defined $2 ? $2 : 1; $maxL = $d if $d > $maxL }
  # Sqrt[3] power
  my $a = ev_part($parts[$k],$E1,$E2,QQ(5),QQ(1)); my $b = ev_part($parts[$k],$E1,$E2,QQ(5),QQ(3));
  my $q = $b/$a; my $p = ($q == 3) ? 1 : ($q == QQ(1)/3 ? -1 : 'X');
  die "unexpected Sqrt3 power" if $p eq 'X';
  my $fac = ($p == 1) ? QQ(1) : QQ(1)/3;   # Sqrt3^p = Sqrt3 * fac
  # Laurent interpolation in s with exponents lo..hi, lo = m - maxL, hi = m (m = largest denominator exponent)
  my $mden = 0; $mden = $_ > $mden ? $_ : $mden for @m;
  my ($lo,$hi) = ($mden - $maxL, $mden); my $n = $hi-$lo+1;
  my @ss = map { QQ($_+2)/QQ(7) } 0..($n+7);
  my @val = map { my $sv = $_; my $Lv = 1/$sv - $al; ev_part($parts[$k],$E1,$E2,$Lv,QQ(1)) * $fac * $sv**(-$lo) } @ss;
  # solve Vandermonde on first n points (Newton), verify on the rest
  my @xs = @ss[0..$n-1]; my @c = @val[0..$n-1];
  for my $j (1..$n-1) { for (my $i=$n-1; $i>=$j; $i--) { $c[$i] = ($c[$i]-$c[$i-1])/($xs[$i]-$xs[$i-$j]); } }
  my @poly = (QQ(0)) x $n; my @basis = (QQ(1));
  for my $i (0..$n-1) { for my $mm (0..$#basis) { $poly[$mm] += $c[$i]*$basis[$mm]; }
    my @nb = (QQ(0)) x (scalar(@basis)+1); for my $mm (0..$#basis) { $nb[$mm+1] += $basis[$mm]; $nb[$mm] -= $xs[$i]*$basis[$mm]; } @basis = @nb; }
  my $ok = 1; for my $r ($n..$#ss) { my $v = QQ(0); $v += $poly[$_]*$ss[$r]**$_ for 0..$n-1; $ok = 0 if $v != $val[$r]; }
  my %S; for my $i (0..$n-1) { $S{$i+$lo} = $poly[$i] if $poly[$i] != 0; }
  my $mn = lp_minexp(\%S); my $mx = lp_maxexp(\%S);
  my $form = ($mn >= 1 && $mx <= $k+1) ? 'FORM OK' : 'FORM FAILS';
  # deg N_k = k+1 - min s-exponent
  my $degN = $k+1-$mn;
  print "k=$k: Sqrt3^$p; denominator (E+3L)^m with m = @m, max L-power elsewhere $maxL; interp verified on 8 extra points: ",
        ($ok?'yes':'NO'), "; s-exponents $mn..$mx; deg N_k = $degN (<= $k: ", ($degN <= $k ? 'yes':'NO'), "); $form\n";
  print "   S_$k(s) = ", lp_str(\%S), "\n";
  print $out "$k ", join(' ', map { "$_:$S{$_}" } sort { $::a <=> $::b } keys %S), "\n";
}
close $out;
