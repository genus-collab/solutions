# test_lib.pl -- sanity checks of hodge_lib.pl against classical values
use strict; use warnings;
use FindBin; require "$FindBin::Bin/hodge_lib.pl";
sub show { my ($name, $got, $want) = @_; printf "%-40s %-20s %-20s %s\n", $name, "$got", "$want", ($got == Math::BigRat->new($want) ? "ok" : "MISMATCH"); }
show("<tau_1>_1", kw(1, 1), "1/24");
show("<tau_4>_2", kw(2, 4), "1/1152");
show("<tau_2 tau_3>_2", kw(2, 2, 3), "29/5760");
show("<tau_7>_3", kw(3, 7), "1/82944");
show("<tau_10>_4", kw(4, 10), "1/7962624");
show("<lambda_1>_{1,1}", hi(1, 1, 0, 0), "1/24");
show("b_2 = <psi^2 lambda_2>_{2,1}", vertex(2, lambda_i(2, 2), [Q(0), Q(0), Q(1)]), "7/5760");
show("b_3 = <psi^4 lambda_3>_{3,1}", vertex(3, lambda_i(3, 3), [Q(0), Q(0), Q(0), Q(0), Q(1)]), "31/967680");
show("Int_M2 lambda_1^3", hi(2, 3, 0), "1/2880");
show("Int_M2 lambda_1 lambda_2", vertex(2, cls_mul(lambda_i(1,2), lambda_i(2,2), 2)), "1/5760");
show("Int_M3 lambda_2^3", vertex(3, cls_mul(cls_mul(lambda_i(2,3), lambda_i(2,3), 3), lambda_i(2,3), 3)), "1/725760");
show("Int_M3 lambda_1^6", hi(3, 6, 0), "1/90720");
# lambda_g lambda_{g-1} formula (FP (9)): Int_{M_{3,1}} psi^2 lambda_3 lambda_2 = 1/2903040 ; Int_{M_{2,1}} psi lambda_2 lambda_1 = 1/2880?
show("Int_M{2,1} psi lambda_2 lambda_1", vertex(2, cls_mul(lambda_i(1,2), lambda_i(2,2), 2), [Q(0), Q(1)]), "1/2880");
show("Int_M{3,1} psi^2 lambda_3 lambda_2", vertex(3, cls_mul(lambda_i(2,3), lambda_i(3,3), 3), [Q(0), Q(0), Q(1)]), "1/120960");
# Mumford check: c(E)c(E^vee) = 1 => lambda_g^2 = 0 ; Int_{M_{2,1}} lambda_2^2 psi^0? (degree 4 on dim 4)
show("Int_M{2,1} lambda_2^2", vertex(2, cls_mul(lambda_i(2,2), lambda_i(2,2), 2), [Q(1)]), "0");
print "Int_M3 lambda_3 lambda_1^3 (no reference value, printed only): ", vertex(3, cls_mul(lambda_i(3,3), cls_mul(lambda_i(1,3), cls_mul(lambda_i(1,3), lambda_i(1,3), 3), 3), 3)), "
";
