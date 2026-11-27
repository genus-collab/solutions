# closed_ref.pl (referee 2): the closed formulas of the TR2 report against the independent recursion RefTR.pm, mod p.
#  (i)   Norlund identities: exp-coefficients of sum a_k t^(2k), sum b_k t^(2k) (first principles) = closed forms of (5.2);
#  (ii)  one-point (5.1)/(5.2) vs TR, g = 1..G, j = 1..2g-2+JX;
#  (iii) (6.1) (via P2 normalisation) vs TR;
#  (iv)  two-point (5.4) vs TR, both orders of (j1,j2), g = 0..G;
#  (v)   general-n integrated formula: n = 2 (consistency with (iv)) and n = 3 vs TR (g <= G-1, small mu).
# usage: perl closed_ref.pl tn td G [JX] [J2max] [J3max] [prime]
use strict; use warnings; use integer; use FindBin; use lib $FindBin::Bin; use RefTR; use RefXYc;
my ($tn,$td,$G,$JX,$J2,$J3,$p) = @ARGV; $JX //= 4; $J2 //= 6; $J3 //= 0; $p //= 2147483647;
my $t0 = time;
my $S = RefTR->new($tn,$td,3*$G-1,$p); $S->run($G,2);
my $tau = $S->{tau}; my $sig = $S->{sig}; my $X = RefXYc->new($tau,$p,$G+1);
my ($nN,$bN,$n1,$b1,$n61,$b61,$n2,$b2,$n3,$b3,$nsym) = (0)x11;
# (i)
for my $j (1..3*$G+3) { my ($A,$B) = $X->AB($j); my ($An,$Bn) = $X->AB_norlund($j);
  for my $m (0..$G+1) { $nN += 2; $bN++ if $A->[$m] != $An->[$m]; $bN++ if $B->[$m] != $Bn->[$m]; } }
# (ii), (iii)
for my $g (1..$G) { for my $j (1..2*$g-2+$JX) { my $tr = $S->xcoef($g,$j); my $f = $X->onept($g,$j); $n1++;
    if ($f != $tr) { $b1++; print "MISMATCH one-point g=$g j=$j\n"; }
    # (6.1): sigma^g Xi_g(j) = prod_{l=1}^{2g-2}(j-l) [lam^(2g)] S1(tau lam)^(1+M) S1((1+tau) lam)^(1-N)
    my $N = $j*$tau % $p; my $M = $j*(($tau+1)%$p) % $p; my $sa = $X->s1pow(RefTR::md(1+$M,$p),$g); my $sb = $X->s1pow(RefTR::md(1-$N,$p),$g);
    my $br = 0; for my $q (0..$g) { $br = ($br + RefTR::pw($tau,2*$q,$p) * $sa->[$q] % $p * RefTR::pw(($tau+1)%$p,2*$g-2*$q,$p) % $p * $sb->[$g-$q]) % $p; }
    my $pre = 1; $pre = $pre * RefTR::md($j - $_,$p) % $p for 1..2*$g-2; my $sx = $pre * $br % $p;    # sigma^g Xi_g(j)
    my $bin = RefTR::gbin(RefTR::md($M-1,$p),$j-1,$p); next unless $bin;
    my $f61 = $sx * RefTR::pw($sig,1-2*$g,$p) % $p * $bin % $p; $n61++; if ($f61 != $tr) { $b61++; print "MISMATCH (6.1) g=$g j=$j\n"; } } }
# (iv)
for my $g (1..$G) { for my $j1 (1..$J2) { for my $j2 ($j1..$J2) { my $tr = $S->xcoef($g,$j1,$j2);
      my $fa = $X->twopt($g,$j1,$j2); my $fb = $X->twopt($g,$j2,$j1); $n2++;
      if ($fa != $tr || $fb != $tr) { $b2++; print "MISMATCH two-point g=$g ($j1,$j2): ".($fa==$tr?"":"order12 ").($fb==$tr?"":"order21")."\n"; }
      if ($g <= 3 && $j2 <= 4) { my $fn = $X->npt($g,$j1,$j2); $nsym++; print "MISMATCH npt(n=2) g=$g ($j1,$j2)\n" if $fn != $tr; $b2++ if $fn != $tr; } } } }
# (v) n = 3
if ($J3) { for my $g (0..$G-1) { for my $a (1..$J3) { for my $b ($a..$J3) { for my $c ($b..$J3) {
          my $tr = $S->xcoef($g,$a,$b,$c); my $f = $X->npt($g,$a,$b,$c); $n3++; if ($f != $tr) { $b3++; print "MISMATCH three-point g=$g ($a,$b,$c)\n"; } } } } } }
print "tau=$tn/$td G=$G p=$p: Norlund $nN checks, $bN mismatches; one-point $n1 checks, $b1 mismatches; (6.1) $n61 checks, $b61 mismatches; ".
  "two-point $n2 checks (+$nsym n=2 general), $b2 mismatches; three-point $n3 checks, $b3 mismatches; ".(time-$t0)."s\n";
