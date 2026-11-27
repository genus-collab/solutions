use strict; use warnings; use FindBin; use lib $FindBin::Bin; use lib "$FindBin::Bin/.."; use Math::BigRat; use SerS;
sub Q { Math::BigRat->new(shift) }
my @bp = (SerS::sone()); for my $p (1..3) { push @bp, SerS::smul($bp[-1], SerS::beta()); }
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
for my $j (0..3) { for my $k (0..2) { my $r = [Q(0)]; for my $p (0..$k) { $r = SerS::padd($r, SerS::pscale(SerS::coeff($bp[$p], $k), Q($j) ** $p / fact($p))); } print "author B^$j [L^$k] = @$r\n"; } }
