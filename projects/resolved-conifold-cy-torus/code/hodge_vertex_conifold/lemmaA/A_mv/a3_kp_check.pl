#!/usr/bin/perl
# a3_kp_check.pl -- (i) the dispersionless Hirota equation e^{D(z)D(w)F_0} = 1 - (d1 D(z)F_0 - d1 D(w)F_0)/(z-w), D(z) = sum_k z^-k d/dP_k,
#   for the genus-0 Marino-Vafa potential F_0(P) = sum_mu sigma^{l-1} Pi_mu |mu|^{l-3} P_mu/|Aut mu|, checked at P = 0 up to weight W;
# (ii) the first KP equation at the dilaton point, with ALL terms (report A3.5):
#   sigma^2 Xi_{G-1}(R+1^4) + 6 sigma^2 sum_{g1+g2=G} sum_{R1+R2=R} Xi_{g1}(R1+1^2) Xi_{g2}(R2+1^2) + 12 sigma Pi_2^2 Xi_G(R+2+2) - 12 sigma Pi_3 Xi_G(R+3+1) = 0.
# Usage: perl a3_kp_check.pl tau W "G:R;G:R;..."
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use mvlib;
my ($tau, $W, $list) = @ARGV; $tau = R($tau); my $sig = $tau*(1+$tau);
sub Pi { my $j = shift; my $p = R(1); for my $b (1..$j-1) { $p *= ($j*$tau + $b); } return $p / fact($j-1); }
# (i) genus 0 at P = 0: second derivatives F_ab = sigma Pi_a Pi_b/(a+b);  LHS = exp(sum F_ab z^-a w^-b), RHS = 1 + sum_k sigma Pi_k/(k+1) sum_{a+b=k+1} z^-a w^-b
my %L; # coefficients of z^-a w^-b of exp(S), S = sum F_ab z^-a w^-b, truncated a,b <= W
my %S; for my $a (1..$W) { for my $b (1..$W) { $S{"$a,$b"} = $sig * Pi($a) * Pi($b) / ($a + $b); } }
my %E = ("0,0" => R(1)); my %Pw = ("0,0" => R(1));
for my $n (1..$W) { my %N; for my $k1 (keys %Pw) { my ($a1,$b1) = split /,/, $k1; for my $k2 (keys %S) { my ($a2,$b2) = split /,/, $k2; my ($a,$b) = ($a1+$a2, $b1+$b2); next if $a > $W || $b > $W; $N{"$a,$b"} = ($N{"$a,$b"} // R(0)) + $Pw{$k1} * $S{$k2}; } }
  %Pw = %N; for my $k (keys %Pw) { $E{$k} = ($E{$k} // R(0)) + $Pw{$k} / fact($n); } }
my $bad = 0; my $cnt = 0;
for my $a (1..$W) { for my $b (1..$W) { next if $a + $b > $W + 1; my $lhs = $E{"$a,$b"} // R(0); my $k = $a + $b - 1; my $rhs = $sig * Pi($k) / ($k + 1); $cnt++; $bad++ if $lhs != $rhs; } }
print "(i) dispersionless Hirota, genus-0 MV at P=0, a+b <= ", $W+1, ": $cnt coefficients, $bad mismatches\n";
# (ii)
my %memo;
sub XiPure { my ($g, @R) = @_; @R = sort { $b <=> $a } @R; my $key = "$g|@R"; return $memo{$key} if exists $memo{$key};
  my $r = scalar @R; my $J = 0; $J += $_ for @R; my $N = 3*$g - 3 + $r; my $M = mv_unit(\@R, $tau, $g, $N);
  my $xi = R(0); for my $a (0..$N) { $xi += R((-1)**$a) / fact($a) * lag($N-$a, 2*$g-2+$r+$a, -$J-$a) * $M->[$g][$a]; }
  $memo{$key} = $xi; return $xi; }
sub Xi { my ($g, @R) = @_; my $J = 0; $J += $_ for @R; my $r = scalar @R; die "empty" if $r == 0;
  if ($g == 0) { my $v = R(1); if ($r >= 3) { $v *= ($J - $_) for 1..($r-3); } elsif ($r == 2) { $v = R(1)/$J; } else { $v = R(1)/($J*($J+1)); } return $v; }
  my @Rp = grep { $_ >= 2 } @R; my $b = $r - scalar @Rp; my $Jp = 0; $Jp += $_ for @Rp;
  return XiPure($g, @Rp) if $b == 0;
  if (@Rp) { my $f = R(1); $f *= ($Jp + $_) for 0..$b-1; return $f * XiPure($g, @Rp); }
  return fact($b-1) * (-(1 + 1/$sig)/24) if $g == 1; return R(0); }
sub subsets { my @a = @_; my @out; for my $m (0..(2**@a - 1)) { my (@I, @J); for my $i (0..$#a) { if (($m >> $i) & 1) { push @I, $a[$i]; } else { push @J, $a[$i]; } } push @out, [\@I, \@J]; } return @out; }
for my $item (split /;/, ($list // '')) { my ($G, $Rs) = split /:/, $item; my @R = grep { $_ ne '' } split /,/, ($Rs // '');
  my $t1 = $G >= 1 ? $sig**2 * Xi($G-1, @R, 1,1,1,1) : R(0);
  my $t2 = R(0); for my $g1 (0..$G) { for my $IJ (subsets(@R)) { $t2 += 6 * $sig**2 * Xi($g1, @{$IJ->[0]}, 1, 1) * Xi($G - $g1, @{$IJ->[1]}, 1, 1); } }
  my $t3 = 12 * $sig * Pi(2)**2 * Xi($G, @R, 2, 2); my $t4 = -12 * $sig * Pi(3) * Xi($G, @R, 3, 1);
  printf "(ii) G=%d R=(%s): total = %s   [terms: %s | %s | %s | %s]\n", $G, join(',', @R), $t1+$t2+$t3+$t4, $t1, $t2, $t3, $t4; }
# (iii) second KP equation at the dilaton point:  sigma Pi_2 Xi_{G-1}(R+2+1^3) + 6 sigma Pi_2 sum Xi_{g1}(R1+1^2) Xi_{g2}(R2+2+1)
#        + 6 Pi_2 Pi_3 Xi_G(R+2+3) - 6 Pi_4 Xi_G(R+4+1) = 0
for my $item (split /;/, ($ARGV[3] // '')) { my ($G, $Rs) = split /:/, $item; my @R = grep { $_ ne '' } split /,/, ($Rs // '');
  my $t1 = $sig * Pi(2) * Xi($G-1, @R, 2, 1,1,1);
  my $t2 = R(0); for my $g1 (0..$G) { for my $IJ (subsets(@R)) { $t2 += 6 * $sig * Pi(2) * Xi($g1, @{$IJ->[0]}, 1, 1) * Xi($G - $g1, @{$IJ->[1]}, 2, 1); } }
  my $t3 = 6 * Pi(2) * Pi(3) * Xi($G, @R, 2, 3); my $t4 = -6 * Pi(4) * Xi($G, @R, 4, 1);
  printf "(iii) G=%d R=(%s): total = %s   [terms: %s | %s | %s | %s]\n", $G, join(',', @R), $t1+$t2+$t3+$t4, $t1, $t2, $t3, $t4; }
