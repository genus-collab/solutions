# ===================== residue.pl (appended to the row construction of edges_sing.pl, colour 0 only) =====================
# Residue lemma test. sigma_a(psi) := Phat_{a,0}(-psi) (rows of the singular colour, frame B, normalised, c1 = 1 and c1 = C1cf).
# Each s-coefficient of each component is fitted EXACTLY to  c0 + Sum_{j=1}^K b_j lam_j(psi),  lam_j(psi) = 1/(E - j psi),
# from z^1..z^K and verified on z^(K+1)..z^NZ.  Then, with C_il := [lam_i lam_l'] P, P(psi,psi') := -<sigma(psi), sigma(psi')>
# (index 0 = the constant function 1), the algebraic identity (proved in the report)
#   [P(psi,psi') - P(-psi',psi')]/(psi+psi') = Sum_{i,l>=1} il/(E(i+l)) C_il lam_i lam_l'  +  Sum_{i>=1} (i/E) S_i lam_i(psi) lam_{-i}(psi'),
#   S_i := E C_i0 + Sum_{l>=1} i C_il/(i+l),   lam_{-i}(psi') = 1/(E + i psi'),
# is compared with the code's E_00 (difference quotient, as edges_sing.pl). Symplecticity (N(psi,-psi) = 0) <=> S_i = 0 (all i) and C_00 = -kappa.
my $NE = $NZ-1;
sub realq { my ($x)=@_; die "non-real coefficient at colour 0" unless $x->[1] == 0; return $x->[0]; }
sub solve { my ($A)=@_; my $n = scalar(@$A); my @M = map { [ @$_ ] } @$A;
  for my $c (0..$n-1) { my $p = $c; $p++ while $p < $n && $M[$p][$c] == 0; die "singular" if $p == $n; @M[$c,$p] = @M[$p,$c];
    my $iv = $ONE/$M[$c][$c]; $M[$c] = [ map { $_*$iv } @{$M[$c]} ];
    for my $r (0..$n-1) { next if $r == $c; my $f = $M[$r][$c]; next if $f == 0; $M[$r] = [ map { $M[$r][$_] - $f*$M[$c][$_] } 0..$n ]; } }
  return [ map { $M[$_][$n] } 0..$n-1 ]; }
sub leaffit { my ($f)=@_;
  for my $K (0..$NZ-1) { my @b = ();
    if ($K > 0) { my @A = map { my $m = $_; [ (map { (QQ(-$_))**$m / $E**($m+1) } 1..$K), $f->[$m] ] } 1..$K; @b = @{ solve(\@A) }; }
    my $c0 = $f->[0]; $c0 -= $b[$_-1]/$E for 1..$K;
    my $ok = 1; for my $m ($K+1..$NZ) { my $v = $ZERO; $v += $b[$_-1]*(QQ(-$_))**$m/$E**($m+1) for 1..$K; $ok = 0 if $v != $f->[$m]; }
    return ($c0, \@b, $K, $NZ-$K) if $ok; }
  return (); }
my %eta = ('0,0'=>$et00, '1,1'=>$et11, '1,2'=>$et12, '2,1'=>$et12);
# code edge E_00 (as edges_sing.pl section (C))
my @Ed00; { my @pre; for my $p (0..$NZ) { for my $q (0..$NZ-$p) { $pre[$p][$q] = sscal(K((($p+$q)%2)?$ONE:-$ONE), pairing($rows[0],$rows[0],$p,$q)); } }
  for my $m (0..$NE) { for my $n (0..$NE-$m) { my $acc = sconst($K0); for my $k (0..$n) { my $x = $pre[$m+1+$k][$n-$k]; $acc = sadd($acc, ($k%2)? sscal(K(-$ONE),$x) : $x); } $Ed00[$m][$n] = $acc; } } }
my $edN = 99; my $edv = 99; for my $m (0..$NE) { for my $n (0..$NE-$m) { my $g = sstrip($Ed00[$m][$n]); $edN = $g->{N} if $g->{N} < $edN; $edv = $g->{v} if $g->{v} < $edv && $g->{v} < $g->{N}; } }
print "code E_00: lowest s-power $edv, known to s^", $edN-1, "\n";
my $nfail = 0;
for my $set (['c1=1', $rows[0]], ['c1=C1cf', $rowsTrue[0]]) {
  my ($name, $R) = @$set; print "==== rows $name ====\n";
  my (%coef, @va, @ea);
  for my $a (0..2) { my $v = 99; my $N = 99;
    for my $m (0..$NZ) { my $g = sstrip($R->[$a][$m]); $N = $g->{N} if $g->{N} < $N; $v = $g->{v} if $g->{v} < $v && $g->{v} < $g->{N}; }
    $va[$a] = $v; $ea[$a] = $N-1;
    for my $e ($v..$N-1) { my @f = map { realq(sget($R->[$a][$_],$e)) } 0..$NZ;
      my ($c0,$b,$Kk,$nch) = leaffit(\@f);
      if (!defined $c0) { print "  sigma_$a [s^$e]: NOT of leaf form c0 + Sum b_j lam_j with K <= ", $NZ-2, "\n"; $nfail++; next; }
      $coef{$a}{$e} = [ $c0, @$b ];
      my $top = 0; for my $j (1..$Kk) { $top = $j if $b->[$j-1] != 0; }
      print "  sigma_$a [s^$e]: c0 = $c0; leaves j <= $top (fit K = $Kk, $nch extra z-orders verified)", ($e >= 0 && $top > $e+1 ? "  [EXCEEDS e+1]" : ""), ($e < 0 && $top > 0 ? "  [psi-DEPENDENT at negative order]" : ""), "\n"; } }
  my $kmin = 99; my $kmax = 99;
  for my $ab (keys %eta) { my ($a,$b) = split /,/, $ab; my $x = $va[$a] + $ea[$b]; $kmax = $x if $x < $kmax; my $y = $va[$a]+$va[$b]; $kmin = $y if $y < $kmin; }
  print "  pairing coefficients complete for s^$kmin .. s^$kmax\n";
  for my $k ($kmin..$kmax) {
    my %C; my $maxi = 0;
    for my $ab (keys %eta) { my ($a,$b) = split /,/, $ab;
      for my $e1 ($va[$a]..$ea[$a]) { my $e2 = $k - $e1; next unless exists $coef{$a}{$e1} && exists $coef{$b}{$e2};
        my ($ca, $cb) = ($coef{$a}{$e1}, $coef{$b}{$e2});
        for my $i (0..$#$ca) { next if $ca->[$i] == 0; for my $l (0..$#$cb) { next if $cb->[$l] == 0;
          $C{$i}{$l} = ($C{$i}{$l} // $ZERO) - $eta{$ab}*$ca->[$i]*$cb->[$l]; $maxi = $i if $i > $maxi; $maxi = $l if $l > $maxi; } } } }
    my $Cg = sub { $C{$_[0]}{$_[1]} // $ZERO };
    # symmetry, constant term, residues S_i
    my $sym = 1; for my $i (0..$maxi) { for my $l (0..$maxi) { $sym = 0 if $Cg->($i,$l) != $Cg->($l,$i); } }
    my $kap = realq(sget($Delta[0],$k));
    my $c00ok = ($Cg->(0,0) == -$kap);
    my @S; my $Sok = 1; for my $i (1..$maxi) { my $v = $E*$Cg->($i,0); $v += QQ($i)*$Cg->($i,$_)/QQ($i+$_) for 1..$maxi; $S[$i] = $v; $Sok = 0 if $v != 0; }
    # predicted edge coefficients A_il and support
    my %A; my @out; for my $i (1..$maxi) { for my $l (1..$maxi) { my $v = QQ($i*$l)/($E*QQ($i+$l))*$Cg->($i,$l); $A{$i}{$l} = $v; push @out, "($i,$l)" if $v != 0 && ($i-1)+($l-1) > $k; } }
    my $anyA = 0; for my $i (1..$maxi) { for my $l (1..$maxi) { $anyA = 1 if $A{$i}{$l} != 0; } }
    # compare with code E_00 at s^k: all (m,n), m+n <= NE; prediction includes the S_i lam_i lam_{-i}' term (zero if symplectic)
    my $cmp = 'not available (edge not known at this order)'; my $ncmp = 0;
    if ($k < $edN) { my $ok = 1;
      for my $m (0..$NE) { for my $n (0..$NE-$m) { my $pr = $ZERO;
        for my $i (1..$maxi) { for my $l (1..$maxi) { next if $A{$i}{$l} == 0; $pr += $A{$i}{$l} * QQ($i)**$m * QQ($l)**$n / $E**($m+$n+2); }
          $pr += (QQ($i)/$E)*$S[$i] * QQ($i)**$m * (QQ(-$i))**$n / $E**($m+$n+2) if $S[$i] != 0; }
        $ncmp++; $ok = 0 if $pr != realq(sget($Ed00[$m][$n],$k)); } }
      $cmp = $ok ? "EQUAL ($ncmp entries m+n <= $NE)" : "DIFFERENT"; $nfail++ unless $ok; }
    $nfail++ unless $sym && $c00ok && $Sok;
    print "  [s^$k] C symmetric: ", ($sym?"yes":"NO"), "; C_00 = -[s^$k]kappa (", $Cg->(0,0), " vs kappa ", $kap, "): ", ($c00ok?"yes":"NO"),
          "; residues S_i = 0 (i=1..$maxi): ", ($Sok?"yes":"NO (".join(',',map { $S[$_] } 1..$maxi).")"), "\n";
    print "         predicted E_00 = Sum il/(E(i+l)) C_il lam_i lam_l': ", ($anyA ? "" : "ZERO; "), "support outside (i-1)+(l-1) <= $k: ", (@out ? join(' ',@out) : "none"), "; vs code E_00: $cmp\n";
    if ($anyA) { my $top = 0; for my $i (1..$maxi) { for my $l (1..$maxi) { $top = $i if $A{$i}{$l} != 0 && $i > $top; } }
      print "         A_il:"; for my $i (1..$top) { print "\n          ", join('  ', map { $A{$i}{$_} } 1..$top); } print "\n"; }
  }
}
print "TOTAL failures: $nfail\ndone\n";
