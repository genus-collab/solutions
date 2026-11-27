# lemmaD_check.pl: direct test of the degree bound (Lemma D) for the exact EO correlators of C_tau, all n:
#   [x^mu] F_{g,n} = 0 for |mu| <= 2g-3+n, and (generically) non-zero at |mu| = 2g-2+n.  Mod two primes (TRmod.pm).
# usage: perl lemmaD_check.pl tn td G NMAX
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use TRmod; use ratrec qw(primes_below);
my ($tn,$td,$G,$NM)=@ARGV; $NM//=3;
my @P = primes_below(2**31,2); my %val;
for my $p (@P){ my $C = TRmod::new_curve($tn,$td,3*$G-3+$NM,$p); TRmod::run_tr($C,$G,$NM);
  for my $g (2..$G){ for my $n (1..$NM){ next unless exists $C->{W}{"$g,$n"}; my $lim=2*$g-2+$n;
      my @tup; my $rec; $rec = sub { my ($pre,$left,$min)=@_; if (@$pre==$n){ push @tup,[@$pre]; return; } for my $m ($min..$left-($n-@$pre-1)){ $rec->([@$pre,$m],$left-$m,$m); } }; $rec->([], $lim, 1);
      for my $t (@tup){ push @{$val{"$g,$n"}{join(",",@$t)}}, TRmod::xcoef($C,$g,@$t); } } } }
my ($tot,$bad)=(0,0);
for my $gn (sort keys %val){ my ($g,$n)=split /,/,$gn; my ($z,$nzl,$viol)=(0,0,0); my @first;
  for my $mk (sort keys %{$val{$gn}}){ my @mu=split /,/,$mk; my $s=0; $s+=$_ for @mu; my $iszero = !grep { $_ } @{$val{$gn}{$mk}}; $tot++;
    if ($s <= 2*$g-3+$n) { if ($iszero) { $z++ } else { $viol++; $bad++; print "VIOLATION ($g,$n) mu=($mk)\n"; } }
    else { $nzl++ unless $iszero; push @first, "($mk)".($iszero?"=0":"!=0"); } }
  print "tau=$tn/$td (g,n)=($g,$n): $z zeros in the range |mu|<=".(2*$g-3+$n).", $viol violations; first layer |mu|=".(2*$g-2+$n).": ".join(" ",@first)."\n"; }
print "tau=$tn/$td: $tot coefficients, $bad violations\n";
