use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
our ($h); BEGIN { $h = shift // 2; }
do "$FindBin::Bin/candK1_lib.pl";
