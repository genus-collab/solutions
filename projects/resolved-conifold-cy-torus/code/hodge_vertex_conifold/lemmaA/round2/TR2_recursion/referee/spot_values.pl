# Referee spot check of the author's exact values (report s.7) with the referee's TR engine and Route E's Hodge tables (mod 2 primes).
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use RefMod; use RefTR; use RefHodge;
my @cases = ( [3,1,5,9,"-558409145/13759414272"], [3,7,5,9,"7833182869/1343692800000"], [-5,3,5,9,"-3452589/204800000"],
              [2,1,4,7,"-182345/23887872"], [2,1,3,5,"-245/82944"], [3,7,2,3,"-3053/2592000"] );
for my $p (prime_list(2)) { setp($p);
  for my $c (@cases) { my ($tn,$td,$g,$j,$val) = @$c; my ($n,$d) = split m{/}, $val;
    my $C = RefTR::new_curve($tn,$td,6*$g-2); RefTR::run($C,$g,1);
    my $tr = RefTR::xi_p2($C,$g,$j); my $hx = RefHodge::xi($g,$tn,$td,$j);
    my $nm = Math::BigInt->new($n)->bmod($p)->numify; my $dm = Math::BigInt->new($d)->bmod($p)->numify; my $v = ratm($nm,$dm);
    print "p=$p tau=$tn/$td Xi_$g($j): author=$val  refTR ".($tr == $v ? "=" : "!=")."  Hodge ".($hx == $v ? "=" : "!=")."\n"; } }
