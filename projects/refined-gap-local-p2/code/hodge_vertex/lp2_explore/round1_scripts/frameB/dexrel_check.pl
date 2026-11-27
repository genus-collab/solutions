use strict; use warnings; use Math::BigRat;
# dexrel_check.pl -- exact y-series check (evidence only) of BS eq:DEXrel (lem:modrel, second relation E = (D-X)/2), whose proof
# in BS v2 (main.tex 2697-2750) is a sketch with details omitted. After cancelling the log y terms, eq:DEXrel reads
#   (2/e1) theta I_1^[1] - f - 9 theta g - A log(1+27y) = 0,
# with f = 3 Sum c_d (-y)^d, c_d = (3d-1)!/(d!)^3, g = Sum c_d (H_{3d-1} - H_d) (-y)^d, A = 1 + theta f,
# (2/e1) theta I_1^[1] = -6 Sum d c_d H_d (-y)^d.  Also checks A = 2F1(1/3,2/3;1;-27y).
my $N = $ARGV[0] // 40;
sub Q { Math::BigRat->new("$_[0]") }
my @H=(Q(0)); for my $n (1..3*$N){ push @H, $H[-1]+Q(1)/Q($n); }
my @fact=(Q(1)); for my $n (1..3*$N){ push @fact, $fact[-1]*Q($n); }
my @c=(Q(0)); for my $d (1..$N){ push @c, $fact[3*$d-1]/($fact[$d]**3); }
my @sg = map { ($_ % 2) ? Q(-1) : Q(1) } 0..$N;
my @A=(Q(1)); my @lhs=(Q(0));
for my $d (1..$N){ push @A, Q(3*$d)*$c[$d]*$sg[$d];
  my $t1 = Q(-6)*Q($d)*$c[$d]*$H[$d]*$sg[$d];
  my $f = Q(3)*$c[$d]*$sg[$d];
  my $tg = Q(9)*Q($d)*$c[$d]*($H[3*$d-1]-$H[$d])*$sg[$d];
  push @lhs, $t1 - $f - $tg; }
# log(1+27y) = Sum_{k>=1} (-1)^(k+1) 27^k y^k / k
my @lg=(Q(0)); for my $k (1..$N){ push @lg, (($k%2)?Q(1):Q(-1))*(Q(27)**$k)/Q($k); }
for my $n (1..$N){ my $s=Q(0); for my $k (1..$n){ $s += $A[$n-$k]*$lg[$k]; } $lhs[$n] -= $s; }
my $bad = grep { $lhs[$_] != 0 } 0..$N;
print "eq:DEXrel (log-free part) vanishes through y^$N: ", ($bad ? "FAIL ($bad coefficients)" : "OK"), "\n";
my @h=(Q(1)); for my $n (1..$N){ push @h, $h[-1]*(Q($n-1)+Q('1/3'))*(Q($n-1)+Q('2/3'))/(Q($n)**2)*Q(-27); }
my $bad2 = grep { $A[$_] != $h[$_] } 0..$N;
print "A = 1 + theta f equals 2F1(1/3,2/3;1;-27y) through y^$N: ", ($bad2 ? "FAIL" : "OK"), "\n";
