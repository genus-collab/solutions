#!/usr/bin/perl
# a3_relation.pl -- check of the "2-insertion relation" obtained from the Marino-Vafa cut-and-join on the unit-leg line at the dilaton
# point s = -1 (report section A3.3):
#   E_g(R) := Xi_g(R+2) - Sum_i gam_{j_i} Xi_g(R_i^+) + Sum_{i<k} gam_{j_i j_k} Xi_g(R_ik)
#             + (sigma/(2 sigma')) Sum_i Sum_{a+b=j_i} a b (Pi_a Pi_b/Pi_{j_i}) [ Xi_{g-1}(a,b,R-j_i) + Sum_{g1+g2=g} Sum_{I+J=R-j_i} Xi_{g1}(a,R_I) Xi_{g2}(b,R_J) ]  = 0,
#   gam_j = (j+1) Pi_{j+1}/(Pi_j sigma'),  gam_{jk} = (j+k) Pi_{j+k}/(Pi_j Pi_k sigma'),  sigma' = 2 tau + 1,  Pi_j = prod_{a=1}^{j-1}(j tau + a)/(j-1)!.
# Conventions: Xi_0(R) = Gamma(J)/Gamma(J-r+3); unit parts by the string equation (Xi_1(1) = -(1+1/sigma)/24 and Xi_g(1^b) = 0, g >= 2).
# Usage: perl a3_relation.pl g tau "R1;R2;..."
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($g0, $tau, $Rlist) = @ARGV; $tau = R($tau);
my @Rs = map { [ grep { $_ ne '' } split /,/, $_ ] } split /;/, $Rlist;
my $sig = $tau*(1+$tau); my $sp = 2*$tau + 1;
sub Pi { my $j = shift; my $p = R(1); for my $b (1..$j-1) { $p *= ($j*$tau + $b); } return $p / fact($j-1); }
my %memo;
sub XiPure { my ($g, @R) = @_; @R = sort { $b <=> $a } @R; my $key = "$g|@R"; return $memo{$key} if exists $memo{$key};
  my $r = scalar @R; my $J = 0; $J += $_ for @R; my $N = 3*$g - 3 + $r; my $M = mv_unit(\@R, $tau, $g, $N);
  my $xi = R(0); for my $a (0..$N) { $xi += R((-1)**$a) / fact($a) * lag($N-$a, 2*$g-2+$r+$a, -$J-$a) * $M->[$g][$a]; }
  $memo{$key} = $xi; return $xi; }
sub Xi { my ($g, @R) = @_; my $J = 0; $J += $_ for @R; my $r = scalar @R; return undef if $r == 0;
  if ($g == 0) { my $v = R(1); if ($r >= 3) { $v *= ($J - $_) for 1..($r-3); } elsif ($r == 2) { $v = R(1)/$J; } else { $v = R(1)/($J*($J+1)); } return $v; }
  my @Rp = grep { $_ >= 2 } @R; my $b = $r - scalar @Rp; my $Jp = 0; $Jp += $_ for @Rp;
  if ($b == 0) { return XiPure($g, @Rp); }
  if (@Rp) { my $f = R(1); $f *= ($Jp + $_) for 0..$b-1; return $f * XiPure($g, @Rp); }
  if ($g == 1) { return fact($b-1) * (-(1 + 1/$sig)/24); }
  return R(0); }
sub subsets { my @a = @_; my @out; for my $m (0..(2**@a - 1)) { my (@I, @J); for my $i (0..$#a) { if (($m >> $i) & 1) { push @I, $a[$i]; } else { push @J, $a[$i]; } } push @out, [\@I, \@J]; } return @out; }
for my $Rref (@Rs) { my @R = @$Rref; my $r = scalar @R; my $g = $g0;
  my $E = Xi($g, @R, 2); my @terms = ("Xi_g(R+2)=" . $E);
  for my $i (0..$r-1) { my $j = $R[$i]; my @Rp = @R; $Rp[$i] = $j + 1; my $gam = ($j+1) * Pi($j+1) / (Pi($j) * $sp); my $v = Xi($g, @Rp); $E -= $gam * $v; push @terms, "-gam_$j*Xi(@Rp)=" . (-$gam*$v); }
  for my $i (0..$r-1) { for my $k ($i+1..$r-1) { my @Rm = map { $R[$_] } grep { $_ != $i && $_ != $k } 0..$r-1; my $jj = $R[$i] + $R[$k];
      my $gam = $jj * Pi($jj) / (Pi($R[$i]) * Pi($R[$k]) * $sp); my $v = Xi($g, @Rm, $jj); $E += $gam * $v; push @terms, "+merge"; } }
  for my $i (0..$r-1) { my $j = $R[$i]; my @rest = map { $R[$_] } grep { $_ != $i } 0..$r-1;
    for my $al (1..$j-1) { my $be = $j - $al; my $cf = $sig / (2*$sp) * $al * $be * Pi($al) * Pi($be) / Pi($j);
      my $t = R(0);
      if ($g - 1 >= 0) { my $v = Xi($g-1, $al, $be, @rest); $t += $v if defined $v; }
      for my $g1 (0..$g) { my $g2 = $g - $g1; for my $IJ (subsets(@rest)) { my ($I, $Jx) = @$IJ; my $v1 = Xi($g1, $al, @$I); my $v2 = Xi($g2, $be, @$Jx); $t += $v1 * $v2; } }
      $E += $cf * $t; } }
  printf "g=%d tau=%s R=(%s): E = %s   [%s]\n", $g, $tau, join(',', @R), $E, join(' ', @terms[0..($#terms < 3 ? $#terms : 3)]); }
