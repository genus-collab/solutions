# tau1_identity.pl (referee 2): claim 9 of the TR2 report at tau = 1,
#    x dF_{g,1}/dx = t Brk_g(v,v) / ((1 - t^2) dX^2)   at v = 1/t,
# with Brk_g(v,v) = omega_{g-1,2}(v,v) + sum_{h=1}^{g-1} omega_{h,1}(v) omega_{g-h,1}(v); tested as power series in t (mod p)
# with the independent recursion RefTR.pm.  (x dF/dx = omega_{g,1}/dX.)
use strict; use warnings; use integer; use FindBin; use lib $FindBin::Bin; use RefTR;
my ($G,$T) = @ARGV; $G //= 6; $T //= 30; my $p = 2147483647;
my $S = RefTR->new(1,1,3*$G-2,$p); $S->run($G,1);
# phi_k(t): xi_k = phi_k(t) dt at v = 1/t:  phi_k = - sum_j c_{k,j} t^j (1 - t)^(-j-2)   (tau = 1)
my @phi; for my $k (0..3*$G) { my @f = (0) x ($T+1);
  for my $j (0..2*$k) { my $c = $S->{c}[$k][$j] or next; for my $i (0..$T-$j) { $f[$j+$i] = RefTR::md($f[$j+$i] - $c * RefTR::gbin($j+1+$i,$i,$p), $p); } } $phi[$k] = \@f; }
# t/q(t) with q = (1 - t)/(1 + t):  t (1+t)/(1-t)
my @toq = (0) x ($T+1); for my $i (1..$T) { $toq[$i] = ($i == 1) ? 1 : 2; }
sub mul { RefTR::smul($_[0],$_[1],$T,$p) }
sub one { my ($g)=@_; my @s = (0) x ($T+1); while (my ($k,$c) = each %{$S->{W}{"$g,1"}}) { $s[$_] = ($s[$_] + $c * $phi[$k][$_]) % $p for 0..$T; } mul(\@s,\@toq) }   # omega_{g,1}/dX
my $bad = 0;
for my $g (2..$G) { my $lhs = one($g);
  my @b = (0) x ($T+1);
  while (my ($key,$c) = each %{$S->{W}{($g-1).",2"}}) { my ($a,$bb) = split /,/, $key; my $pr = mul($phi[$a],$phi[$bb]); $b[$_] = ($b[$_] + $c * $pr->[$_]) % $p for 0..$T; }
  my $brk = mul(mul(\@b,\@toq),\@toq);                      # omega_{g-1,2}(v,v)/dX^2
  for my $h (1..$g-1) { my $pr = mul(one($h),one($g-$h)); $brk->[$_] = ($brk->[$_] + $pr->[$_]) % $p for 0..$T; }
  # rhs = t/(1-t^2) * brk
  my @f = (0) x ($T+1); for my $i (0..$T) { $f[$i] = 1 if $i % 2 == 1; } my $rhs = mul($brk,\@f);
  my @d = grep { $lhs->[$_] != $rhs->[$_] } 0..$T; $bad += @d;
  my ($lo) = grep { $lhs->[$_] } 0..$T;
  print "tau=1 g=$g: identity ".(@d ? "FAILS at t^(@d)" : "holds through t^$T").", lowest power of x dF/dx: t^".($lo // "none")." (Lemma D: >= ".(2*$g-1).")\n"; }
print "tau=1: claim-9 identity, g = 2..$G: ".($bad ? "$bad failures" : "holds")."\n";
