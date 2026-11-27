# closed_n34.pl (referee 2): the general-n integrated formula of report s.5 (iterated residues over n-cycles) for n = 3, 4
# against the independent recursion RefTR.pm, mod p.
# usage: perl closed_n34.pl tn td G3 J3 G4 J4 [prime]
use strict; use warnings; use integer; use FindBin; use lib $FindBin::Bin; use RefTR; use RefXYc;
my ($tn,$td,$G3,$J3,$G4,$J4,$p) = @ARGV; $p //= 2147483647; my $t0 = time;
my $Gr = $G3 >= $G4+1 ? $G3 : $G4+1; my $G = $Gr; my $S = RefTR->new($tn,$td,3*$Gr+1,$p); $S->run($Gr, 3);
my $X = RefXYc->new($S->{tau},$p,$G+2);
my ($n3,$b3,$n4,$b4,$nz3,$nz4) = (0)x6;
for my $g (0..$G3) { for my $a (1..$J3) { for my $b ($a..$J3) { for my $c ($b..$J3) {
        my $tr = $S->xcoef($g,$a,$b,$c); my $f = $X->npt($g,$a,$b,$c); $n3++; $nz3++ if $tr;
        if ($f != $tr) { $b3++; print "MISMATCH n=3 g=$g ($a,$b,$c)\n"; } } } } }
for my $g (0..$G4) { for my $a (1..$J4) { for my $b ($a..$J4) { for my $c ($b..$J4) { for my $d ($c..$J4) {
          my $tr = $S->xcoef($g,$a,$b,$c,$d); my $f = $X->npt($g,$a,$b,$c,$d); $n4++; $nz4++ if $tr;
          if ($f != $tr) { $b4++; print "MISMATCH n=4 g=$g ($a,$b,$c,$d)\n"; } } } } } }
print "tau=$tn/$td: n=3 (g<=$G3, mu_i<=$J3) $n3 checks ($nz3 non-zero), $b3 mismatches; n=4 (g<=$G4, mu_i<=$J4) $n4 checks ($nz4 non-zero), $b4 mismatches; ".(time-$t0)."s\n";
