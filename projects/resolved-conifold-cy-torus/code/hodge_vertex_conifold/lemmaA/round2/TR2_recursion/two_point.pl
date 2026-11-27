# two_point.pl: closed two-point function of C_tau from ABDKS (17), n = 2, after integrating by parts in both variables:
#   [x1^j1 x2^j2] F_{g,2} = (1/(j1 j2)) [hbar^(2g)] Res_{z1=oo} Res_{z2=oo} Psi_{j1}(z1) Psi_{j2}(z2) dz1 dz2 / ((z1-z2)^2 - c^2),
#   c = (j1+j2) hbar/2, inner residue in z2 taken with |z2| > |z1| (expansion of (z2-z1)^(-k) in z1/z2),
#   Psi_j(z) = exp(-j S(j hbar d) xhat) = w^(-j) (1+w)^(j(tau+1)) Sum_m (hbar w)^(2m) Sum_p A_{m-p}(j) B_p(j) (1+w)^(-2p),  w = 1/z.
# Compared with exact TR (TRmod over 3 primes, rational reconstruction).  usage: perl two_point.pl tn td G JMAX
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use TRmod; use ratrec qw(primes_below recon);
sub R { Math::BigRat->new(@_) }
my ($tn,$td,$G,$JM)=@ARGV; $JM//=2*$G+1; my $tau=R("$tn/$td");
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
sub binomr { my ($a,$i)=@_; return R(0) if $i<0; my $r=R(1); for my $k (0..$i-1){ $r *= ($a-$k); } $r/fact($i) }
my $K=$G+2; my @s = map { R(1)/(R(4)**$_ * fact(2*$_+1)) } 0..$K; my @sp=(R(1)); for my $n (1..$K){ my $acc=R(0); $acc += $s[$_]*$sp[$n-$_] for 1..$n; push @sp, -$acc; }
sub ck { my ($u,$a,$k)=@_; my $c=R(0); for my $p (0..$k){ $c += $s[$p]*$sp[$k-$p]*($u**(2*$p))*($a**(2*($k-$p))); } $c }
my @al=(R(1)/$tau, R(-1)/($tau+1));
sub expco { my ($c,$m)=@_; my @ex=(R(1)); for my $q (1..$m){ my $acc=R(0); $acc += $_*$c->[$_]*$ex[$q-$_] for 1..$q; push @ex, $acc/$q; } \@ex }
my %AB; sub AB { my ($j,$G)=@_; return @{$AB{$j}} if $AB{$j}; my @Ak=(R(0)); my @Bk=(R(0));
  for my $k (1..$G){ push @Ak, fact(2*$k-1)*($j/$al[0])*ck(R($j),$al[0],$k); push @Bk, fact(2*$k-1)*($j/$al[1])*ck(R($j),$al[1],$k); }
  $AB{$j}=[expco(\@Ak,$G),expco(\@Bk,$G)]; @{$AB{$j}} }
# coefficient [w^e] of (1+w)^M w^(2m) P_m(w),  P_m = Sum_p A_{m-p} B_p (1+w)^(-2p)
sub psicoef { my ($j,$m,$e)=@_; my ($A,$B)=AB($j,$G); my $M=$j*($tau+1); my $q=$e-2*$m; return R(0) if $q<0; my $s=R(0);
  for my $p (0..$m){ $s += $A->[$m-$p]*$B->[$p]*binomr($M-2*$p,$q); } $s }
sub two_closed { my ($g,$j1,$j2)=@_; my $tot=R(0); my $cc=R($j1+$j2)/2;
  # [hbar^(2g)]: hbar^(2m1) from Psi1, hbar^(2m2) from Psi2, hbar^(2mk) from the kernel, m1+m2+mk = g
  for my $mk (0..$g){ for my $m1 (0..$g-$mk){ my $m2=$g-$mk-$m1; my $k=2*$mk+2;
      for my $l (0..$j2){ my $e2 = $j2-1-2*$mk-$l; last if $e2<0; # Psi2 coefficient [w2^(j2+1-k-l)] = [w2^e2]
        my $c1 = psicoef($j1,$m1,1+$l+$j1); next if $c1->is_zero; my $c2 = psicoef($j2,$m2,$e2); next if $c2->is_zero;
        $tot += ($cc**(2*$mk))*binomr(R($k+$l-1),$l)*$c1*$c2; } } }
  $tot/($j1*$j2) }
my @P = primes_below(2**31,3); my %val;
for my $p (@P){ my $C=TRmod::new_curve($tn,$td,3*$G-1,$p); TRmod::run_tr($C,$G,2);
  for my $g (1..$G){ for my $j1 (1..$JM){ for my $j2 ($j1..$JM){ push @{$val{"$g,$j1,$j2"}}, TRmod::xcoef($C,$g,$j1,$j2); } } } }
my ($n,$bad)=(0,0);
for my $g (1..$G){ my @row; for my $j1 (1..$JM){ for my $j2 ($j1..$JM){ my $f=two_closed($g,$j1,$j2); my $r=recon($val{"$g,$j1,$j2"},\@P); $n++;
      my $okp=1; for my $i (0..$#P){ my $p=$P[$i]; my $fm = TRmod::ratmod($f->numerator->copy->bmod($p)->numify, $f->denominator->copy->bmod($p)->numify, $p); $okp=0 unless $fm == $val{"$g,$j1,$j2"}[$i]; }
      # exact agreement mod each of the 3 primes (rational reconstruction alone may fail for large rationals)
      my $ok = $okp; $bad++ unless $ok; push @row, "($j1,$j2)".($f->is_zero?"=0":"").($ok?"":"[MISMATCH f=$f tr=".(defined $r?$r:"?")."]"); } }
  print "tau=$tn/$td g=$g: ".join(" ",@row)."\n"; }
print "tau=$tn/$td: $n two-point coefficients compared, $bad mismatches\n";
