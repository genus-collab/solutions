# status_check61.pl -- evaluate the conjectural one-leaf formula (6.1) of proof_M_residue.md,
#   sigma^h Xi_h(j) = prod_{l=1}^{2h-2} (j-l) * [lambda^{2h}] S(e1 l)^{1+j e2} S(e2 l)^{1+j e1},
#   e1 = -tau, e2 = 1+tau, sigma = tau(1+tau), S(x) = (x/2)/sinh(x/2),
# and print the prediction in the TR normalisation of review_checks/referee1_tr*.pl:  TR = sigma^(1-h) Xi_h  (c = 1).
# Usage: perl status_check61.pl TAU H JMAX [odd]
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
my ($tau, $h, $jmax, $odd) = (R(shift // 1), shift // 5, shift // 13, shift // 0);
my @B = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"));  # B_{2k}, k = 0..7
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
my $K = $h;
my @lS = (R(0)); for my $k (1..$K){ push @lS, -$B[$k]/(2*$k*fact(2*$k)); }
sub pexp { my @c=@_; my @r=(R(1)); push @r, R(0) for 1..$K;
  for my $n (1..$K){ my $s=R(0); for my $m (1..$n){ $s += $m*$c[$m]*$r[$n-$m]; } $r[$n]=$s/$n; } return @r; }
my $e1 = -$tau; my $e2 = 1+$tau; my $sig = $tau*(1+$tau);
for my $jj (1..$jmax) { next if $odd && $jj % 2 == 0; my $j = R($jj);
  my @c=(R(0)); for my $k (1..$K){ push @c, (1+$j*$e2)*$lS[$k]*($e1**(2*$k)) + (1+$j*$e1)*$lS[$k]*($e2**(2*$k)); }
  my @r = pexp(@c); my $br = $r[$h];
  my $p = R(1); $p *= ($j - $_) for 1..2*$h-2;
  my $xi = $p*$br/($sig**$h); my $tr = $xi*($sig**(1-$h));
  print "tau=$tau h=$h j=$jj: bracket=", $br->bstr, "  Xi(c=1)=", $xi->bstr, "  TR-normalised=", $tr->bstr, "\n"; }
