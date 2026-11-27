my $c1sec = sscal(K(-$ONE), $wc);   # C1cf = -wc
# ---- h_w^-1 (frame B) ----
my @B = (QQ('1/6'),QQ('-1/30'),QQ('1/42'),QQ('-1/30'),QQ('5/66'));
my @w = ($E, -$e1, -$e2); my @lg = ($ZERO) x ($NZ+1);
for my $k (1..5) { my $m = 2*$k-1; next if $m > $NZ; my $S = $ZERO; $S += $ONE/($_**$m) for @w; $lg[$m] = $B[$k-1]/(2*$k*(2*$k-1))*$S; }
my @hinv = ($ONE, ($ZERO) x $NZ);
for my $n (1..$NZ) { my $v = $ZERO; for my $k (1..$n) { $v += $k*$lg[$k]*$hinv[$n-$k]; } $hinv[$n] = $v/$n; }
# ---- rows per colour ----
my @om = (K($ONE), $W, $W2);
my (@rows, @Delta, @rowsTrue);
for my $i (0..2) {
  my $oi = $om[$i]; my $ti = spoly(0, [ ksc(QQ('-1/3'),$oi), ksc(QQ('1/3'),$oi) ]);        # t_i = omega^i (s-1)/3 (exact)
  my ($vi, $ivi);
  if ($i == 0) { $vi = spoly(1, [ K($ONE) ]); $ivi = spoly(-1, [ K($ONE) ]); }            # v_0 = 1+3t_0 = s exactly
  else { $vi = spoly(0, [ ksub(K($ONE),$oi), $oi ]); $ivi = sinv($vi); }
  my $Li = sscal(K(-$alpha), $ivi);
  my $it = sinv($ti);
  # (alpha+L)^-1 = v/(3 alpha t): exact in terms of v and 1/t (no Laurent cancellation needed)
  my $iApL = sscal(K($ONE/(3*$alpha)), smul($vi, $it));
  my $tov = smul($ti, $ivi);
  my $Dop = sub { my ($F)=@_; sadd(sadd(smul(sscal(K(QQ('1/3'))
            , sadd($sS, sconst(K(-$ONE)))), sder($F)), sscal(K(QQ('-1/3')),$F)), smul($tov,$F)) };
  # L^j (alpha+L)^-k = (-alpha)^j v^(k-j) (3 alpha t)^-k : every term regular, no precision loss
  my @Vp = (sconst(K($ONE))); push @Vp, smul($Vp[-1],$vi) for 1..$NZ;
  my @Tp = (sconst(K($ONE))); push @Tp, smul($Tp[-1],$it) for 1..$NZ;
  my @P0; for my $k (0..$NZ) { my $pp = (3+$pk[$k])/2 - $k; my $pref = $pp >= 0 ? QQ(3)**$pp : $ONE/(QQ(3)**(-$pp)); $pref *= $alpha*$alpha/((3*$alpha)**$k);
    my $acc = sconst($K0); for my $j (0..$k) { $acc = sadd($acc, sscal(K($Nc[$k][$j]*(-$alpha)**$j), $Vp[$k-$j])); }
    push @P0, sscal(K($pref), smul($acc, $Tp[$k])); }
  my @P1; for my $m (0..$NZ) { my $x = smul($Li,$P0[$m]); $x = sadd($x, $Dop->($P0[$m-1])) if $m >= 1; push @P1, $x; }
  my @P2; for my $m (0..$NZ) { my $x = smul($Li,$P1[$m]); $x = sadd($x, ssub($Dop->($P1[$m-1]), smul($XX,$P1[$m-1]))) if $m >= 1; push @P2, smul($Del,$x); }
  # frame B
  my @R; for my $a ([@P0],[@P1],[@P2]) { my @r; for my $k (0..$NZ) { my $acc = sconst($K0); for my $j (0..$k) { $acc = sadd($acc, sscal(K($hinv[$k-$j]), $a->[$j])); } push @r, $acc; } push @R, \@r; }
  $rows[$i] = \@R;
  $Delta[$i] = sscal(K(3*$e1*$e2*$E**3), smul(smul($ti,$ti), smul($ivi,$ivi)));
  my $ic1 = sinv($c1sec);
  $rowsTrue[$i] = [ $R[0], [ map { smul($R[1][$_], $ic1) } 0..$NZ ],
     [ map { sadd(smul($c1sec,$R[2][$_]), smul(sscal(K($E/2), ssub($c1sec,$ic1)), $R[1][$_])) } 0..$NZ ] ];
}
my $et00 = 3*$alpha**3*$e1*$e2; my $et11 = -3*$E*$e1*$e2; my $et12 = -3*$e1*$e2;
