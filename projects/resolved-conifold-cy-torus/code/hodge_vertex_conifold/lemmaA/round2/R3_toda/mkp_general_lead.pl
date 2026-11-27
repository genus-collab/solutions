#!/usr/bin/perl
# Leading order at the dilaton point of the u_k-linear mKP equation P_k(D) tau_0.tau_1 = 0 (mkp_eqs.pm), for every
# target profile R (parts >= 2) and x-power n:  Sum_terms c Sum_{assignments of the labelled parts of R to the factors}
#   Sum_{(g_a,k_a): Sum_a (2 g_a - 2 + m_a + k_a) = n} prod_a f_a = 0,
# f_a = [2 if S-factor and k_a = 0] sigma^{g-1+m} tau^k (prod_{i in I} Pi_i) L(g,k,U,u)/k!,  U = T_a + I_a, m = |I_a| + u_a,
# L = J_U (J_U+1)...(J_U+u-1) (delta_g(U))_k X_g(U)  (U nonempty),  delta_g(U) = D(U) - 2g + 2,
# L(empty): g = 0: k = 0 -> (u-1)!/2, k = 1,2 -> (u-1)!, k >= 3 -> 0;  g = 1: k = 0 -> (u-1)! X_1(1), else 0;  g >= 2: 0.
# Usage: perl mkp_general_lead.pl k Gmin Gmax rmax Jmax [taus]   (PERTURB="g|parts|delta" as in mkp_lead.pl)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use XiData; use mkp_eqs;
$| = 1;
sub R { Math::BigRat->new(@_) }
my ($K, $Gmin, $Gmax, $rmax, $Jmax, $taustr) = @ARGV; $taustr //= "2/7,3/2,-5/3";
my @taus = map { R($_) } split /,/, $taustr;
my $terms = mkp_eqs::eq_k($K);
print "equation k=$K: ", mkp_eqs::show($terms), "\n";
my ($PG, $PT, $PD) = split /\|/, ($ENV{PERTURB} // "");
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub ff { my ($x, $k) = @_; my $p = R(1); for my $i (0..$k-1) { $p *= ($x - $i); } $p }
sub rf { my ($x, $k) = @_; my $p = R(1); for my $i (0..$k-1) { $p *= ($x + $i); } $p }
sub Pi { my ($k, $tau) = @_; my $p = R(1); for my $a (1..$k-1) { $p *= ($k*$tau + $a); } $p / fact($k-1) }
sub X { my ($g, $T, $tau) = @_; my $v = XiData::xi_at($g, $T, $tau);
  if (defined $PD && $g == $PG && join(',', sort { $a <=> $b } @$T) eq $PT) { $v += R($PD); } $v }
sub L { my ($g, $k, $U, $u, $tau) = @_; my $sig = $tau*($tau+1);
  if (!@$U) { return R(0) if $u < 1;
    if ($g == 0) { return fact($u-1)/2 if $k == 0; return fact($u-1) if $k <= 2; return R(0); }
    if ($g == 1) { return $k == 0 ? fact($u-1) * (-(1 + 1/$sig)/24) : R(0); }
    return R(0); }
  my $J = 0; $J += $_ for @$U; my $D = $J - scalar(@$U); my $f = ff(R($D - 2*$g + 2), $k); return R(0) if $f == 0;
  return rf(R($J), $u) * $f * X($g, $U, $tau); }
sub fval { my ($fac, $g, $k, $T, $tau) = @_; my ($type, $ins, $u) = @$fac; my $sig = $tau*($tau+1); my $m = scalar(@$ins) + $u;
  return R(0) if $type eq 'A' && $k < 1;
  my $Lv = L($g, $k, [@$T, @$ins], $u, $tau); return R(0) if $Lv == 0;
  my $c = ($sig ** ($g - 1 + $m)) * ($tau ** $k) / fact($k); $c *= 2 if $type eq 'S' && $k == 0;
  for my $i (@$ins) { $c *= Pi($i, $tau); }
  return $c * $Lv; }
sub relation { my ($n, $Rr, $tau) = @_; my @Rr = @$Rr; my $r = scalar @Rr; my $tot = R(0);
  for my $t (@$terms) { my ($c, $fac) = @$t; my $q = scalar @$fac;
    my $nass = $q ** $r;
    for my $code (0..$nass-1) { my @T = map { [] } 1..$q; my $cc = $code; for my $i (0..$r-1) { push @{$T[$cc % $q]}, $Rr[$i]; $cc = int($cc / $q); }
      # enumerate (g_a, k_a)
      my $rec; $rec = sub { my ($a, $budget, $acc) = @_;   # budget = remaining Sum (2g + k - 2 + m)
        my $fa = $fac->[$a]; my $m = scalar(@{$fa->[1]}) + $fa->[2]; my $kmin = $fa->[0] eq 'A' ? 1 : 0;
        if ($a == $q - 1) { my $rem = $budget + 2 - $m;   # = 2g + k
          for my $g (0..int($rem/2)) { my $k = $rem - 2*$g; next if $k < $kmin; my $v = fval($fa, $g, $k, $T[$a], $tau); next if $v == 0; $tot += $c * $acc * $v; }
          return; }
        my $minrest = 0; for my $b ($a+1..$q-1) { my $fb = $fac->[$b]; $minrest += scalar(@{$fb->[1]}) + $fb->[2] + ($fb->[0] eq 'A' ? 1 : 0) - 2; }
        my $maxcont = $budget - $minrest;   # cont_a = 2g + k - 2 + m <= maxcont
        for my $g (0..int(($maxcont + 2 - $m)/2)) { for my $k ($kmin..($maxcont + 2 - $m - 2*$g)) { my $cont = 2*$g + $k - 2 + $m;
            my $v = fval($fa, $g, $k, $T[$a], $tau); next if $v == 0; $rec->($a + 1, $budget - $cont, $acc * $v); } } };
      $rec->(0, $n, R(1)); } }
  return $tot; }
my @Rs = ([]); my $rr; $rr = sub { my ($pre, $min, $J) = @_; for (my $j = $min; $J + $j <= $Jmax; $j++) { my @p = (@$pre, $j); push @Rs, [@p]; $rr->(\@p, $j, $J + $j) if scalar(@p) < $rmax; } };
$rr->([], 2, 0);
my ($ntot, $nbad) = (0, 0);
for my $G ($Gmin..$Gmax) { for my $Rr (@Rs) { for my $par (0, 1) { my $n = 2*$G + $par;
  my @RK = (@$Rr, $K); my $DK = 0; $DK += $_ - 1 for @RK;
  my @res = map { relation($n, $Rr, $_) } @taus; my $ok = !(grep { $_ != 0 } @res); $ntot++; $nbad++ unless $ok;
  printf "k=%d G=%d R={%s} %s [R+%d: D=%d excess=%d] %s\n", $K, $G, join(',', @$Rr), ($par ? "ODD " : "EVEN"), $K, $DK, 2*$G-2-$DK, ($ok ? "OK" : "FAIL @res"); } } }
printf "SUMMARY k=%d: %d relations, %d failures\n", $K, $ntot, $nbad;
