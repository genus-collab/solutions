use strict; use warnings; use FindBin; use lib "$FindBin::Bin/.."; use XY qw(xi_formula_n1); use Math::BigRat;
my $X = xi_formula_n1(Math::BigRat->new(1),6,12);
for my $gj ("5,8","5,9","6,10","6,11"){ my ($g,$j)=split /,/,$gj; my $v=$X->{$gj}; my $trn = $v/(Math::BigRat->new(2)**($g-1)); print "tau=1 Xi_$g($j) = $v  (TR normalisation sigma^(1-g) x this = $trn)\n"; }
