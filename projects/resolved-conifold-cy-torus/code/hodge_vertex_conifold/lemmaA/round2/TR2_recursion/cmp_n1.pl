# cmp_n1.pl: compare the x-y swap closed formula (XY.pm, exact BigRat) with exact TR on C_tau (TRmod.pm over several primes,
# CRT + rational reconstruction) for the one-point functions Xi_g(mu), 2 <= g <= G, 1 <= mu <= M_g = 2g-2+JX.
# usage: perl cmp_n1.pl tn td G [JX] [nprimes]
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use TRmod; use ratrec qw(primes_below recon); use XY qw(xi_formula_n1); use Math::BigRat;
my ($tn,$td,$G,$JX,$np)=@ARGV; $JX//=4; $np//=4;
my $tau = Math::BigRat->new("$tn/$td"); my $t0=time;
my $Mmax = 2*$G-2+$JX;
my $XF = xi_formula_n1($tau,$G,$Mmax); my $tf=time-$t0;
my @P = primes_below(2**31,$np); my %vals;
for my $p (@P){ my $C = TRmod::new_curve($tn,$td,3*$G-2,$p); TRmod::run_tr($C,$G,1);
  for my $g (2..$G){ for my $mu (1..2*$g-2+$JX){ push @{$vals{"$g,$mu"}}, TRmod::xi_tr($C,$g,$mu); } } }
my $ttr=time-$t0-$tf; my ($nc,$nbad,$nzero)=(0,0,0);
for my $g (2..$G){ my @row;
  for my $mu (1..2*$g-2+$JX){ my $r = recon($vals{"$g,$mu"},\@P); my $f = $XF->{"$g,$mu"}; $nc++;
    # check the formula value mod each prime as well (independent of reconstruction)
    my $okp = 1; for my $i (0..$#P){ my $p=$P[$i]; my $fm = TRmod::ratmod($f->numerator->copy->bmod($p)->numify, $f->denominator->copy->bmod($p)->numify, $p); $okp=0 unless $fm == $vals{"$g,$mu"}[$i]; }
    my $same = (defined $r && $r == $f && $okp) ? 1 : 0; $nbad++ unless $same; $nzero++ if $f->is_zero;
    push @row, "$mu:".($f->is_zero ? "0" : "nz").($same ? "" : "[MISMATCH tr=".(defined $r ? $r : "?")." f=$f]"); }
  print "tau=$tn/$td g=$g: ", join(" ",@row), "\n"; }
print "tau=$tn/$td: $nc comparisons, $nbad mismatches; formula ${tf}s, TR ${ttr}s; first nonzero values: ";
for my $g (2..$G){ print "Xi_$g(".(2*$g-1).")=".$XF->{"$g,".(2*$g-1)}."  "; } print "\n";
