# Frobenius series at Y=1 of theta^3 - Y prod_{i=1..3}(theta + a_i), a_i = -lam/(3z) + (i-1)/3, theta = Y d/dY, w = Y-1.
use strict; use warnings;
my $lam = 1.7;   # leaf scale c (generic)
sub act { # apply operator to w^rho, return hash offset->coef (power rho+offset)
  my ($rho,$z)=@_; my @a = map { -$lam/(3*$z) + $_/3 } (0,1,2);
  my $th = sub { my ($f)=@_; my %g; for my $k (keys %$f){ my $e=$rho+$k; $g{$k-1}+=$e*$f->{$k}; $g{$k}+=$e*$f->{$k}; } \%g };
  my $f={0=>1}; my $t3 = $th->($th->($th->($f)));
  my $p = {0=>1}; for my $ai (@a){ my $q=$th->($p); $q->{$_}+= $ai*$p->{$_} for keys %$p; $p=$q; }
  my %r = %$t3; for my $k (keys %$p){ $r{$k}-=$p->{$k}; $r{$k+1}-=$p->{$k}; }
  return \%r; }
sub C { my ($rho,$z,$j)=@_; my $r=act($rho,$z); return $r->{$j-2}//0; } # coefficient of w^{rho-2+j}
my $z=0.37; my $s=1+$lam/$z;
printf "indicial check: C0(rho)/(rho(rho-1)(rho-s)) at rho=0.3,2.2: %.10f %.10f\n", map { C($_,$z,0)/($_*($_-1)*($_-$s)) } (0.3,2.2);
for my $j (0..5){ printf "C_%d(0.3) = %.6g\n",$j,C(0.3,$z,$j); }
sub series { my ($z,$N,$rho0)=@_; my @b=(1); for my $n (1..$N){ my $acc=0; for my $j (1..$n){ $acc += C($rho0+$n-$j,$z,$j)*$b[$n-$j]; } $b[$n] = -$acc/C($rho0+$n,$z,0); } @b }
# pole order near z0 = -lam/j for the exponent-s solution
if(!$ENV{SKIP}){for my $j (1..4){ my $z0=-$lam/$j; for my $n (1..4){ my @o; for my $d (1e-3,1e-4){ my $zz=$z0+$d; my @b=series($zz,$n,1+$lam/$zz); push @o, abs($b[$n]); } my $ord = log($o[1]/$o[0])/log(10); printf "j=%d n=%d |b_n| growth exponent ~ %.2f\n",$j,$n,$ord; } }}
if($ENV{SMALLZ}){ for my $n (1..5){ my @o; for my $zz (1e-3,1e-4){ my @b=series($zz,$n,1+$lam/$zz); push @o, abs($b[$n]); } printf "n=%d small-z exponent of |b_n| ~ z^-%.2f\n",$n, log($o[1]/$o[0])/log(10); } }
