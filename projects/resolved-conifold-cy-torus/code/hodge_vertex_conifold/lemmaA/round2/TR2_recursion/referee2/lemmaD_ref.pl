# lemmaD_ref.pl (referee 2): direct test of Lemma D / Theorem 1 with the independent engine RefTR.pm.
# For every (g,n) with g >= 2, n <= NM + G - g computed by one run, and every sorted mu (mu_i >= 1):
#   |mu| <= 2g-3+n  ->  [x^mu] F_{g,n} must vanish (tested mod each prime);   |mu| = 2g-2+n (first layer): reported.
# usage: perl lemmaD_ref.pl tn td G NM [nprimes]
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefTR;
my ($tn,$td,$G,$NM,$np) = @ARGV; $np //= 2;
my @P = (2147483647, 2147483629, 2147483587)[0..$np-1];
my (%val,%first); my $t0 = time;
for my $p (@P) { my $S = RefTR->new($tn,$td,3*$G-3+$NM,$p); $S->run($G,$NM);
  for my $g (2..$G) { for my $n (1..$NM+$G-$g) { my $lim = 2*$g-2+$n; my @tup; my $rec;
      $rec = sub { my ($pre,$left,$min)=@_; if (@$pre == $n) { push @tup, [@$pre]; return; }
        for my $m ($min..$left-($n-@$pre-1)) { last if $m*($n-@$pre) > $left; $rec->([@$pre,$m],$left-$m,$m); } };
      $rec->([], $lim, 1);
      for my $t (@tup) { push @{$val{"$g,$n"}{join(",",@$t)}}, $S->xcoef($g,@$t); } } } }
my ($tot,$bad,$nz1,$n1) = (0,0,0,0);
for my $gn (sort { my @a=split /,/,$a; my @b=split /,/,$b; $a[0]<=>$b[0] || $a[1]<=>$b[1] } keys %val) { my ($g,$n) = split /,/, $gn;
  my ($z,$viol,$fz,$fn) = (0,0,0,0);
  for my $mk (keys %{$val{$gn}}) { my @mu = split /,/, $mk; my $s = 0; $s += $_ for @mu; my $iszero = !grep { $_ } @{$val{$gn}{$mk}};
    if ($s <= 2*$g-3+$n) { $tot++; if ($iszero) { $z++ } else { $viol++; $bad++; print "VIOLATION tau=$tn/$td (g,n)=($g,$n) mu=($mk)\n"; } }
    else { $n1++; if ($iszero) { $fz++ } else { $fn++; $nz1++; } } }
  print "tau=$tn/$td (g,n)=($g,$n): $z in-range coefficients vanish, $viol violations; first layer |mu|=".(2*$g-2+$n).": $fn non-zero, $fz zero\n"; }
print "tau=$tn/$td G=$G NM=$NM primes=$np: $tot in-range coefficients, $bad violations; first layer $nz1/$n1 non-zero; ".(time-$t0)."s\n";
