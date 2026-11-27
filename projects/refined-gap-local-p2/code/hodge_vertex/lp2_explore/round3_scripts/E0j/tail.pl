
# ===================== E0j tail (round 3) =====================
# Checks of the mixed-edge statements against the code's edges (evidence only).
#  (M1) leaf decomposition of the singular column: [s^k] Phat_{a,0}(-psi) = c0 + Sum_{l<=k+1} b_l lam_l(psi), lam_l = 1/(E - l psi)
#       (k = -1: pure constant), in both gauges (c1 = 1 code rows and true c1 = C1cf rows); verified on the spare z-orders.
#  (M2) closed form  E_0j(psi,w) = - Sum_{k1+k2=k} Sum_l l < S_l^(k1), v_{j,k2}(-w) > lam_l(psi) lam_l(-w)   at every order s^k,
#       compared entry by entry with the code's difference quotient Ed{0}{j}[m][n] (constants c0 and the s^-1 part dropped).
#  (M3) transposition E_j0(n,m) = E_0j(m,n) for j = 1,2 and E_21(n,m) = E_12(m,n).
#  (M4) Galois: conj(E_01) = E_02, conj(E_12) = E_21, conj(E_11) = E_22 (conj: omega -> omega^2).
my $KMAX = $ENV{KMAX} // 4;
sub kconj { my ($x)=@_; [ $x->[0]-$x->[1], -$x->[1] ] }
sub solve { # Gaussian elimination: A (n x n of K), rhs (n of K) -> x or undef
  my ($A,$r)=@_; my $n = scalar(@$r); my @M = map { [ @{$A->[$_]}, $r->[$_] ] } 0..$n-1;
  for my $c (0..$n-1) { my $p = $c; $p++ while $p < $n && kz($M[$p][$c]); return undef if $p >= $n; @M[$c,$p] = @M[$p,$c];
    my $iv = kinv($M[$c][$c]); $M[$c] = [ map { kmul($iv,$_) } @{$M[$c]} ];
    for my $q (0..$n-1) { next if $q == $c; my $f = $M[$q][$c]; next if kz($f); $M[$q] = [ map { ksub($M[$q][$_], kmul($f,$M[$c][$_])) } 0..$n ]; } }
  [ map { $M[$_][$n] } 0..$n-1 ] }
my $fails = 0;
my @etadiag; # eta^{ab}
my @ET = ([K($et00),$K0,$K0],[$K0,K($et11),K($et12)],[$K0,K($et12),$K0]);
for my $gauge ('c1=1','c1=C1cf') {
  my $R0 = $gauge eq 'c1=1' ? $rows[0] : $rowsTrue[0];
  # (M1) leaf fits: $S{$k}[a] = { c0 => K, b => [undef, b_1, ..., b_K] }
  my %S;
  for my $k (-1..$KMAX) { my $Kk = $k+1;
    for my $a (0..2) { my @G = map { ksc((($_%2)?-$ONE:$ONE), sget($R0->[$a][$_],$k)) } 0..$NZ;   # psi^m coefficient of Phat_a(-psi)
      my @A; my @r; for my $m (0..$Kk) { my @row = ( $m == 0 ? K($ONE) : $K0 ); push @row, K(QQ($_)**$m/($E**($m+1))) for 1..$Kk; push @A, \@row; push @r, $G[$m]; }
      my $x = solve(\@A,\@r); if (!defined $x) { print "  (M1) $gauge s^$k comp $a: singular system\n"; $fails++; next; }
      my $ok = 1; my $nchk = 0;
      for my $m ($Kk+1..$NZ) { my $v = ($m == 0 ? $x->[0] : $K0); for my $l (1..$Kk) { $v = kadd($v, ksc(QQ($l)**$m/($E**($m+1)), $x->[$l])); } $nchk++; $ok = 0 unless kz(ksub($v,$G[$m])); }
      $S{$k}[$a] = { c0 => $x->[0], b => $x };
      print "  (M1) $gauge [s^$k] Phat_$a(-psi) = c0 + Sum_{l<=$Kk} b_l lam_l: ", ($ok ? "OK" : "FAILED"), " ($nchk spare z-orders); c0 = ", kstr($x->[0]),
        ($Kk ? "; b = (" . join(', ', map { kstr($x->[$_]) } 1..$Kk) . ")" : ""), "\n";
      $fails++ unless $ok; } }
  # (M2) closed form vs code
  for my $j (1,2) { my $Rj = $gauge eq 'c1=1' ? $rows[$j] : $rowsTrue[$j]; my $bad = 0; my $cnt = 0;
    for my $k (0..$KMAX) { last if $k > $NS-3;
      for my $m (0..$NE) { for my $n (0..$NE-$m) { my $pred = $K0;
        for my $k1 (0..$k) { my $k2 = $k - $k1;
          for my $l (1..$k1+1) {
            # <S_l^(k1), v_{j,k2}(-w)>/(E + l w): w^n coefficient
            my $coef = $K0;
            for my $q (0..$n) { my $pv = $K0;
              for my $a (0..2) { for my $b (0..2) { next if kz($ET[$a][$b]); my $sa = $S{$k1}[$a]{b}[$l] // $K0;
                $pv = kadd($pv, kmul(kmul($ET[$a][$b], $sa), ksc((($q%2)?-$ONE:$ONE), sget($Rj->[$b][$q],$k2)))); } }
              my $r = $n - $q; $coef = kadd($coef, ksc(((($r%2)?-$ONE:$ONE) * QQ($l)**$r)/($E**($r+1)), $pv)); }
            $pred = kadd($pred, ksc(-QQ($l) * QQ($l)**$m/($E**($m+1)), $coef)); } }
        $cnt++; $bad++ unless kz(ksub($pred, sget($Ed{0}{$j}[$m][$n],$k))); } } }
    print "  (M2) $gauge: closed form of E_0$j (constants and s^-1 part dropped) vs code difference quotient, s^0..s^", ($KMAX < $NS-3 ? $KMAX : $NS-3), ", m+n <= $NE: ",
      ($bad ? "FAILED ($bad of $cnt)" : "EQUAL ($cnt entries)"), "\n"; $fails += $bad; }
}
# (M2') the s^-1 part of the singular column contributes nothing: [s^-1] N0j(psi,w) is psi-independent
# (M3) transpositions
for my $pr ([1,0],[2,0],[2,1]) { my ($i,$j) = @$pr; my $ok = 1; my $cnt = 0;
  for my $m (0..$NE) { for my $n (0..$NE-$m) { for my $e (0..$NS-3) { $cnt++; $ok = 0 unless kz(ksub(sget($Ed{$i}{$j}[$n][$m],$e), sget($Ed{$j}{$i}[$m][$n],$e))); } } }
  print "  (M3) E_$i$j(n,m) = E_$j$i(m,n), s^0..s^", $NS-3, ": ", ($ok ? "OK" : "NO"), " ($cnt)\n"; $fails++ unless $ok; }
# (M4) Galois conjugation
for my $pr ([[0,1],[0,2]], [[1,0],[2,0]], [[1,2],[2,1]], [[1,1],[2,2]], [[0,0],[0,0]]) { my ($p,$q) = @$pr; my $ok = 1; my $cnt = 0;
  for my $m (0..$NE) { for my $n (0..$NE-$m) { for my $e (0..$NS-3) { $cnt++; $ok = 0 unless kz(ksub(kconj(sget($Ed{$p->[0]}{$p->[1]}[$m][$n],$e)), sget($Ed{$q->[0]}{$q->[1]}[$m][$n],$e))); } } }
  print "  (M4) conj(E_$p->[0]$p->[1]) = E_$q->[0]$q->[1]: ", ($ok ? "OK" : "NO"), " ($cnt)\n"; $fails++ unless $ok; }
# (M5) regular-regular edges: lowest s-power (expect >= 0) -- already in (C); recheck E_12 explicitly and print E_12(0,0) at s^0
print "  (M5) E_12(0,0)|_{s=0} = ", kstr(sget($Ed{1}{2}[0][0],0)), ", E_11(0,0)|_{s=0} = ", kstr(sget($Ed{1}{1}[0][0],0)), "\n";
print "TOTAL failures: $fails\n";
print "done\n";
