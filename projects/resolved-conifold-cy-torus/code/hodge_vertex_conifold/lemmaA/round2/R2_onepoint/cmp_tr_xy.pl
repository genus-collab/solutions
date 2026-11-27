#!/usr/bin/perl
# cmp_tr_xy.pl TRFILE [nmax] : compare the exact EO recursion output of status_trC1.pl (tau = 1; lines "g=G n=N: (mu)=value ...",
# value = [x^mu] F_{g,n} / prod C(2 mu_i - 1, mu_i - 1)) with the x-y swap closed formula (xy_npoint.pl conventions, tau = 1).
use strict; use warnings;
use FindBin; use lib "$FindBin::Bin/lib";
use Math::BigRat; use MP;
my ($file, $nmax) = @ARGV; $nmax //= 4;
my $tau = Q(1); my $e2 = Q(2);
sub fact { my $k = shift; my $r = Q(1); $r *= $_ for 2..$k; return $r; }
sub Pi { my $j = shift; my $p = Q(1); $p *= ($j * $tau + $_) for 1..$j-1; return $p / fact($j-1); }
our ($n, $H, $HT, $NV);
sub qser { my ($mu, $T) = @_; my $n1 = $tau * $mu; my $n2 = $e2 * $mu;
  my @P = map { [ map { Q(0) } 0..$HT ] } 0..$T; $P[0][0] = Q(1);
  my $mulf = sub { my ($F) = @_; my @C = map { [ map { Q(0) } 0..$HT ] } 0..$T;
    for my $i (0..$T) { for my $k (0..$T-$i) { for my $p (0..$HT) { next if $P[$i][$p]->is_zero; for my $q (0..$HT-$p) { next if $F->[$k][$q]->is_zero; $C[$i+$k][$p+$q] += $P[$i][$p] * $F->[$k][$q]; } } } }
    @P = @C; };
  for my $kk (1..$n2) { my $m = Q($kk) - ($n2 + 1)/2; my @F = map { [ map { Q(0) } 0..$HT ] } 0..$T; $F[0][0] = Q(1); if ($T >= 1) { $F[1][0] = Q(1); $F[1][1] = $m / $e2; } $mulf->(\@F); }
  for my $kk (1..$n1) { my $m = Q($kk) - ($n1 + 1)/2; my @F = map { [ map { Q(0) } 0..$HT ] } 0..$T; my $cf = Q(1);
    for my $d (0..$T) { $F[$d][$d] = $cf if $d <= $HT; $cf *= -$m / $tau; } $mulf->(\@F); }
  return \@P; }
sub trunc_h { my ($p) = @_; my %q; for my $k (keys %$p) { my @e = split /,/, $k; next if $e[$n] > $HT; $q{$k} = $p->{$k}; } return \%q; }
sub hpow_poly { my ($arr) = @_; my %p; for my $d (0..$#$arr) { next if $arr->[$d]->is_zero; my @e = (0) x $NV; $e[$n] = $d; $p{join(',', @e)} = $arr->[$d]->copy; } return \%p; }
sub ncycles { my ($k) = @_; return ([0]) if $k == 1; my @out;
  my $rec; $rec = sub { my ($path, $rest) = @_; if (!@$rest) { my @s; for my $i (0..$#$path) { $s[$path->[$i]] = $path->[($i+1) % @$path]; } push @out, \@s; return; }
    for my $i (0..$#$rest) { my @r2 = @$rest; my ($x) = splice(@r2, $i, 1); $rec->([@$path, $x], \@r2); } };
  $rec->([0], [1..$k-1]); return @out; }
sub xyval { my ($g, @mu) = @_; local $n = scalar @mu; local $H = 2*$g - 2 + $n; local $HT = $H; local $NV = $n + 1;
  my $sum = 0; $sum += $_ for @mu; my $T = 2*$sum + 3*$n + 3; my @q = map { qser($_, $T) } @mu; my $tot = Q(0);
  if ($n == 1) { my $mu = $mu[0]; return $q[0][$mu+1][$H+1] / $mu / ($mu) if $H + 1 <= $HT; # need hb^{H+1}
  }
  for my $sg (ncycles($n)) { my @fac; my $sign = 1;
    for my $i (0..$n-1) { my $s = $sg->[$i]; my $k = Q($mu[$i] + $mu[$s]) / 2; if ($i < $s) { push @fac, [$i, $s, $k]; } else { push @fac, [$s, $i, -$k]; $sign = -$sign; } }
    my $P = { join(',', (0) x $NV) => Q($sign) }; my @remaining = @fac;
    for my $step (0..$n-1) { my @here = grep { $_->[0] == $step } @remaining; @remaining = grep { $_->[0] != $step } @remaining; my $m = scalar @here;
      my @A = map { my $f = $_; mp_add(mp_var($f->[1], $NV), mp_scal(-$f->[2], mp_var($n, $NV))) } @here;
      my $emax = 0; for my $key (keys %$P) { my @e = split /,/, $key; $emax = $e[$step] if $e[$step] > $emax; } my $smax = $emax + $mu[$step] + 1;
      my @hs; if ($m == 0) { @hs = (mp_const(1,$NV), map { {} } 1..$smax); } elsif ($m == 1) { @hs = (mp_const(1,$NV)); push @hs, trunc_h(mp_mul($hs[-1], $A[0])) for 1..$smax; }
      else { my @p0 = (mp_const(1,$NV)); my @p1 = (mp_const(1,$NV)); for (1..$smax) { push @p0, trunc_h(mp_mul($p0[-1], $A[0])); push @p1, trunc_h(mp_mul($p1[-1], $A[1])); }
        for my $s (0..$smax) { my $h = {}; $h = mp_add($h, trunc_h(mp_mul($p0[$_], $p1[$s-$_]))) for 0..$s; push @hs, $h; } }
      my $NP = {};
      for my $key (keys %$P) { my @e = split /,/, $key; my $ee = $e[$step]; $e[$step] = 0; my $mono = { join(',', @e) => $P->{$key} };
        for my $a (0..($ee + $mu[$step] + 1 - $m)) { my $s = $ee + $mu[$step] + 1 - $m - $a; next if $s < 0 || $s > $smax || $a > $T;
          my $qa = hpow_poly($q[$step][$a]); next unless %$qa; next unless %{$hs[$s]}; $NP = mp_add($NP, trunc_h(mp_mul(mp_mul($mono, $qa), $hs[$s]))); } }
      $P = $NP; }
    for my $key (keys %$P) { my @e = split /,/, $key; $tot += $P->{$key} if $e[$n] == $H; } }
  my $pm = Q(1); $pm *= $_ for @mu; return ($n % 2 ? 1 : -1) * $tot / $pm; }
open my $fh, '<', $file or die "no $file";
my ($nt, $nb) = (0, 0);
while (my $line = <$fh>) { next unless $line =~ /^g=(\d+) n=(\d+): (.*)$/; my ($g, $nn, $rest) = ($1, $2, $3); next if $nn > $nmax;
  my ($c, $b) = (0, 0);
  while ($rest =~ /\(([\d,]+)\)=(-?[\d\/]+)/g) { my @mu = split /,/, $1; my $v = Q($2);
    my $F;
    if ($nn == 1) { local $HT = 2*$g; local $H = 2*$g - 1; local $n = 1; local $NV = 2; my $mu = $mu[0]; my $q = qser($mu, $mu + 2); $F = $q->[$mu+1][2*$g] / ($mu * $mu); }
    else { $F = xyval($g, @mu); }
    my $pp = Q(1); $pp *= Pi($_) for @mu; my $mine = $F / $pp;
    $c++; if ($mine != $v) { $b++; print "  MISMATCH g=$g mu=(@mu): TR $v  xy $mine\n"; } }
  $nt += $c; $nb += $b; print "g=$g n=$nn: compared $c, mismatches $b\n"; }
print "TOTAL compared $nt, mismatches $nb\n";
