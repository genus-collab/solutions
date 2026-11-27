# t_mod.pl: run TRmod at tau = tn/td over several primes, reconstruct Xi_g(j) (n = 1) rationally.
# usage: perl t_mod.pl tn td G [nprimes] [jextra]
use strict; use warnings; use FindBin; use lib $FindBin::Bin;
use TRmod; use ratrec qw(primes_below recon);
my ($tn,$td,$G,$np,$jx)=@ARGV; $np//=3; $jx//=3;
my @P = primes_below(2**31,$np); my %vals; my $t0=time;
for my $p (@P){ my $C = TRmod::new_curve($tn,$td,3*$G-2,$p); print STDERR "p=$p setup ".(time-$t0)."s\n";
  TRmod::run_tr($C,$G,1,verbose=>($p==$P[0]));
  for my $g (2..$G){ for my $j (1..2*$g-2+$jx){ push @{$vals{"$g,$j"}}, TRmod::xi_tr($C,$g,$j); } } }
print "tau=$tn/$td, primes ".join(",",@P).", time ".(time-$t0)."s\n";
for my $g (2..$G){ my @out; for my $j (1..2*$g-2+$jx){ my $r = recon($vals{"$g,$j"},\@P); push @out, "Xi_$g($j)=".(defined $r ? $r : "?"); } print join("  ",@out),"\n"; }
