package Hodge;
# Exact Hodge integrals on Mbar_{g,n} (no topological recursion; Mumford's GRR formula + Witten-Kontsevich numbers):
#   HK(g, m, beta, a) = Int_{Mbar_{g,n}} ch_1^{m1} ch_3^{m3} ch_5^{m5} kappa_{beta_1}...kappa_{beta_s} psi_1^{a_1}...psi_n^{a_n}
# Mumford:  ch_l(E) = B_{l+1}/(l+1)! [ kappa_l - Sum_i psi_i^l + (1/2) xi_*( Sum_{p+q=l-1} (-1)^p psi_*^p psi_o^q ) ]  (l odd),
# the boundary push-forward running over the irreducible divisor (Mbar_{g-1,n+2}) and all ordered reducible splittings.
# Restrictions to the boundary: ch_l -> ch_l (irreducible) or ch_l' + ch_l'' (reducible); kappa_b likewise (AC kappa classes).
# kappa classes are removed by kappa_b = pi_*(psi_{n+1}^{b+1}), pi^* kappa_b' = kappa_b' - psi_{n+1}^{b'}.
# The CY Hodge class (c = 1, weights 1, tau, -1-tau) is  Omega = Lambda^vee(1) prod_i (1 + v alpha_i(1+alpha_i)),  v = -1/sigma,
# so that Int Omega (...) = Sum_k v^k Xi^(k) directly gives the tau-independent components (proof_status s.1.1(e)).
use strict; use warnings;
use Math::BigRat; use WK;
sub R { Math::BigRat->new(shift) }
our $NCH = $ENV{HODGE_NCH} // 3;   # number of odd Chern characters ch_1, ch_3, ch_5, (ch_7): ch_l with l <= h suffice in genus h (3 through genus 6)
my @CHL = (1, 3, 5, 7, 9)[0..$NCH-1];
my %BERN = (2 => R("1/6"), 4 => R("-1/30"), 6 => R("1/42"), 8 => R("-1/30"), 10 => R("5/66"), 12 => R("-691/2730"));
sub fact { my $n = shift; my $r = R(1); $r *= $_ for 2..$n; return $r; }
my @MCOEF = map { $BERN{$_ + 1} / fact($_ + 1) } @CHL;   # B_{l+1}/(l+1)!
sub binom { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; my $r = 1; for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i } return $r; }

our %memo; our %memo0;
my $ZERO = R(0);
# sub-multisets of a sorted list, returned as [ [sub], [complement], multiplicity ]
sub subms { my @X = @_; my %t; $t{$_}++ for @X; my @vals = sort { $a <=> $b } keys %t; my @out = ([[], [], 1]);
  for my $v (@vals) { my $mlt = $t{$v}; my @new;
    for my $o (@out) { for my $i (0..$mlt) { push @new, [[@{$o->[0]}, ($v) x $i], [@{$o->[1]}, ($v) x ($mlt - $i)], $o->[2] * binom($mlt, $i)]; } }
    @out = @new; }
  return @out; }

sub HK0 {   # Int_{Mbar_{g,n}} kappa_beta psi^a  (no Hodge classes)
  my ($g, $beta, $av) = @_;
  my $n = scalar @$av; return $ZERO if 2*$g - 2 + $n <= 0;
  my $deg = 0; $deg += $_ for @$beta, @$av; return $ZERO if $deg != 3*$g - 3 + $n;
  my @b = sort { $a <=> $b } @$beta; my @aa = sort { $a <=> $b } @$av;
  return WK::wk($g, @aa) unless @b;
  my $key = "$g|" . join(',', @b) . "|" . join(',', @aa);
  return $memo0{$key} if exists $memo0{$key};
  my $b0 = pop @b;   # remove kappa_{b0}
  my $res = R(0);
  for my $s (subms(@b)) { my ($A, $rest, $mult) = @$s; my $sa = 0; $sa += $_ for @$A;
    my $sign = (scalar(@$A) % 2) ? -1 : 1;
    my $v = HK0($g, $rest, [@aa, $b0 + 1 + $sa]); next if $v == 0;
    $res += $v * ($sign * $mult); }
  $memo0{$key} = $res; return $res; }

sub HK {
  my ($g, $m, $beta, $av) = @_;   # $m = [m1, m3, m5]
  my $n = scalar @$av;
  return $ZERO if $g < 0 || 2*$g - 2 + $n <= 0;
  for (@$av) { return $ZERO if $_ < 0 }
  my $dm = 0; for my $i (0..$NCH-1) { $dm += $CHL[$i] * $m->[$i] }
  my $deg = $dm; $deg += $_ for @$beta, @$av;
  return $ZERO if $deg != 3*$g - 3 + $n;
  return HK0($g, $beta, $av) if $dm == 0;
  return $ZERO if $g == 0;
  return $ZERO if $g == 1 && $dm > 1;
  my @b = sort { $a <=> $b } @$beta; my @aa = sort { $a <=> $b } @$av;
  my $key = "$g|" . join(',', @$m) . "|" . join(',', @b) . "|" . join(',', @aa);
  return $memo{$key} if exists $memo{$key};
  my $li; for (my $i = $NCH - 1; $i >= 0; $i--) { if ($m->[$i] > 0) { $li = $i; last } }
  my $l = $CHL[$li]; my @mp = @$m; $mp[$li]--;
  my $tot = R(0);
  # kappa_l term
  $tot += HK($g, \@mp, [@b, $l], \@aa);
  # psi terms
  my %t; $t{$_}++ for @aa;
  for my $d (keys %t) { my @c = @aa; for my $i (0..$#c) { if ($c[$i] == $d) { $c[$i] += $l; last } }
    my $v = HK($g, \@mp, \@b, \@c); $tot -= $v * $t{$d} if $v != 0; }
  # irreducible boundary
  my $half = R("1/2");
  if ($g >= 1) { my $s = R(0);
    for my $p (0..$l-1) { my $v = HK($g - 1, \@mp, \@b, [@aa, $p, $l - 1 - $p]); $s += (($p % 2) ? -$v : $v) if $v != 0; }
    $tot += $half * $s; }
  # reducible boundary (ordered splittings)
  my $s = R(0);
  my @msplits = ([[], [], 1]);
  for my $i (0..$NCH-1) { my @new; for my $o (@msplits) { for my $k (0..$mp[$i]) { push @new, [[@{$o->[0]}, $k], [@{$o->[1]}, $mp[$i] - $k], $o->[2] * binom($mp[$i], $k)]; } } @msplits = @new; }
  my @asub = subms(@aa); my @bsub = subms(@b);
  for my $g1 (0..$g) { my $g2 = $g - $g1;
    for my $as (@asub) { my ($S, $Sc, $am) = @$as; my $n1 = scalar(@$S) + 1; my $n2 = scalar(@$Sc) + 1;
      next if 2*$g1 - 2 + $n1 <= 0 || 2*$g2 - 2 + $n2 <= 0;
      my $sS = 0; $sS += $_ for @$S;
      for my $ms (@msplits) { my ($m1, $m2, $mm) = @$ms;
        my $d1 = 0; for my $i (0..$NCH-1) { $d1 += $CHL[$i] * $m1->[$i] }
        next if $g1 == 0 && $d1 > 0; next if $g2 == 0 && ($dm - $l - $d1) > 0;
        for my $bs (@bsub) { my ($B1, $B2, $bm) = @$bs; my $sB = 0; $sB += $_ for @$B1;
          my $p = 3*$g1 - 3 + $n1 - $d1 - $sB - $sS; next if $p < 0 || $p > $l - 1;
          my $v1 = HK($g1, $m1, $B1, [@$S, $p]); next if $v1 == 0;
          my $v2 = HK($g2, $m2, $B2, [@$Sc, $l - 1 - $p]); next if $v2 == 0;
          my $w = $v1 * $v2 * ($am * $mm * $bm); $s += (($p % 2) ? -$w : $w); } } } }
  $tot += $half * $s;
  my $res = $tot * $MCOEF[$li];
  $memo{$key} = $res; return $res; }

# ---------- the Hodge class in the ch-basis ----------
# polynomials in ch_1, ch_3, ch_5: hash { "m1,m3,m5" => coefficient } ; coefficients may be v-polynomials (array refs)
sub chmul { my ($p, $q, $maxdeg) = @_; my %r;
  for my $x (keys %$p) { my @mx = split /,/, $x; for my $y (keys %$q) { my @my = split /,/, $y;
      my @mz = map { $mx[$_] + $my[$_] } 0..$NCH-1; my $d = 0; $d += $CHL[$_] * $mz[$_] for 0..$NCH-1; next if $d > $maxdeg;
      my $k = join(',', @mz); my $c = vmul($p->{$x}, $q->{$y}); $r{$k} = exists $r{$k} ? vadd($r{$k}, $c) : $c; } }
  return \%r; }
sub vmul { my ($a, $b) = @_; $a = [$a] unless ref $a eq 'ARRAY'; $b = [$b] unless ref $b eq 'ARRAY'; my @c = map { R(0) } 0..($#$a + $#$b);
  for my $i (0..$#$a) { for my $j (0..$#$b) { $c[$i+$j] += $a->[$i] * $b->[$j] } } return \@c; }
sub vadd { my ($a, $b) = @_; $a = [$a] unless ref $a eq 'ARRAY'; $b = [$b] unless ref $b eq 'ARRAY'; my $n = $#$a > $#$b ? $#$a : $#$b;
  return [ map { ($a->[$_] // R(0)) + ($b->[$_] // R(0)) } 0..$n ]; }
sub vscale { my ($a, $c) = @_; $a = [$a] unless ref $a eq 'ARRAY'; return [ map { $_ * $c } @$a ]; }
# lambda_k as ch-polynomials (genus h, degree <= maxdeg): c(E) = exp(Sum_{l odd} (l-1)! ch_l)
sub lambdas { my ($h, $maxdeg) = @_;
  my %X; for my $i (0..$NCH-1) { my $l = $CHL[$i]; next if $l > $maxdeg; my @e = (0) x $NCH; $e[$i] = 1; $X{join(',', @e)} = [fact($l - 1)]; }
  # exp(X) by series
  my %E = (join(',', (0) x $NCH) => [R(1)]); my %term = %E;
  for my $k (1..$maxdeg) { my $t = chmul(\%term, \%X, $maxdeg); %term = map { $_ => vscale($t->{$_}, R(1)/$k) } keys %$t;
    for (keys %term) { $E{$_} = exists $E{$_} ? vadd($E{$_}, $term{$_}) : $term{$_}; } }
  my @lam; for my $k (0..$h) { my %L; for my $x (keys %E) { my @mx = split /,/, $x; my $d = 0; $d += $CHL[$_] * $mx[$_] for 0..$NCH-1; $L{$x} = $E{$x} if $d == $k; } push @lam, \%L; }
  return @lam; }
# Omega = Lambda^vee(1) * Sum_{k,l} lambda_k lambda_l (-r+)^k (-r-)^l,  r+ + r- = -v, r+ r- = v
sub omega_ch { my ($h, $maxdeg) = @_;
  my @lam = lambdas($h, $maxdeg);
  my %LV; for my $k (0..$h) { for (keys %{$lam[$k]}) { my $c = vscale($lam[$k]{$_}, (($k % 2) ? R(-1) : R(1))); $LV{$_} = exists $LV{$_} ? vadd($LV{$_}, $c) : $c; } }
  # Q_n = r+^n + r-^n as v-polynomials: Q_0 = 2, Q_1 = -v, Q_n = -v Q_{n-1} - v Q_{n-2}
  my @Q = ([R(2)], [R(0), R(-1)]); for my $n (2..2*$h) { push @Q, vadd(vmul([R(0), R(-1)], $Q[$n-1]), vmul([R(0), R(-1)], $Q[$n-2])); }
  my %C2;
  for my $k (0..$h) { for my $l ($k..$h) {
      my $coef; # v-poly
      if ($k == $l) { $coef = [ (map { R(0) } 1..$k), R(1) ]; }               # (r+ r-)^k = v^k
      else { my $sgn = (($k + $l) % 2) ? R(-1) : R(1); $coef = vscale(vmul([ (map { R(0) } 1..$k), R(1) ], $Q[$l - $k]), $sgn); }
      my $pr = chmul($lam[$k], $lam[$l], $maxdeg);
      for (keys %$pr) { my $c = vmul($pr->{$_}, $coef); $C2{$_} = exists $C2{$_} ? vadd($C2{$_}, $c) : $c; } } }
  return chmul(\%LV, \%C2, $maxdeg); }

# K = exp(-Sum_k kappa_k/k): hash { "sorted beta" => coefficient } up to degree maxdeg
sub partitions { my ($n, $max) = @_; $max //= $n; return ([]) if $n == 0; my @out;
  for my $k (reverse 1..($max < $n ? $max : $n)) { for my $p (partitions($n - $k, $k)) { push @out, [$k, @$p]; } } return @out; }
sub kexp { my ($maxdeg) = @_; my %K;
  for my $d (0..$maxdeg) { for my $p (partitions($d)) { my %mlt; $mlt{$_}++ for @$p; my $c = R(1);
      for my $k (keys %mlt) { $c *= (R(-1)/$k) ** $mlt{$k} / fact($mlt{$k}); }
      $K{join(',', sort { $a <=> $b } @$p)} = $c; } }
  return \%K; }

# T_h[beta] = Int_{Mbar_h} Omega K kappa_beta  (v-polynomial)
our %Tmemo; our ($OM, $KX, $HCUR); our %Wmemo; our %T1memo;
sub setup { my $h = shift; $HCUR = $h; my $D = 3*$h - 3; $D = 1 if $h == 1; $OM = omega_ch($h, $D); $KX = kexp($h == 1 ? 1 : $D); %Tmemo = (); %Wmemo = (); %T1memo = (); }
# genus 1: the forgetful recursion stops at Mbar_{1,1}; base T1(beta; a) = Int_{Mbar_{1,1}} Omega K kappa_beta psi^a, directly
sub T1 { my ($beta, $ax) = @_; my @beta = sort { $a <=> $b } @$beta; my $key = join(',', @beta) . "|$ax"; return $T1memo{$key} if exists $T1memo{$key};
  my $res = [R(0)]; my $db = 0; $db += $_ for @beta;
  for my $mk (keys %$OM) { my @m = split /,/, $mk;
    for my $gk (keys %$KX) { my @gm = $gk eq '' ? () : split /,/, $gk;
      my $v = HK(1, \@m, [@beta, @gm], [$ax]); next if $v == 0; $res = vadd($res, vscale($OM->{$mk}, $v * $KX->{$gk})); } }
  $T1memo{$key} = $res; return $res; }
sub T { my @beta = sort { $a <=> $b } @_; my $h = $HCUR; my $key = join(',', @beta);
  return $Tmemo{$key} if exists $Tmemo{$key};
  my $db = 0; $db += $_ for @beta; my $D = 3*$h - 3; my $res = [R(0)];
  if ($db <= $D) {
  for my $mk (keys %$OM) { my @m = split /,/, $mk; my $dm = 0; $dm += $CHL[$_] * $m[$_] for 0..$NCH-1; next if $dm + $db > $D;
    for my $gk (keys %$KX) { my @gm = $gk eq '' ? () : split /,/, $gk; my $dg = 0; $dg += $_ for @gm; next if $dm + $db + $dg != $D;
      my $v = HK($h, \@m, [@beta, @gm], []); next if $v == 0;
      $res = vadd($res, vscale($OM->{$mk}, $v * $KX->{$gk})); } } }
  $Tmemo{$key} = $res; return $res; }

# W(beta; a) = Int_{Mbar_{h,n}} Omega K kappa_beta prod psi_i^{a_i}, by forgetting the last point (see report s.2):
# W(n; beta; a) = Sum_{A in beta} [W(n-1; beta-A+{e-1}; a') - W(n-1; beta-A+{e}; a')] + [a_n = 0] Sum_{i<n} W(n-1; beta; a'-e_i),
# e = a_n + |A|, kappa_{-1} = 0, kappa_0 = 2h-2+(n-1).
# (%Wmemo declared above)
sub W { my ($beta, $av) = @_; my $h = $HCUR; my @aa = @$av; my $n = scalar @aa;
  my @b = sort { $a <=> $b } @$beta; @aa = sort { $b <=> $a } @aa;
  return T1(\@b, $aa[0]) if $h == 1 && $n == 1;
  return T(@b) if $n == 0;
  my $key = join(',', @b) . "|" . join(',', @aa);   # a kept in order: last point is forgotten first (symmetric anyway)
  return $Wmemo{$key} if exists $Wmemo{$key};
  my $an = pop @aa; my $res = [R(0)];
  for my $s (subms(@b)) { my ($A, $rest, $mult) = @$s; my $sa = 0; $sa += $_ for @$A; my $e = $an + $sa;
    # + kappa_{e-1}
    if ($e - 1 == 0) { $res = vadd($res, vscale(W($rest, \@aa), R(2*$h - 2 + $n - 1) * $mult)); }
    elsif ($e - 1 > 0) { $res = vadd($res, vscale(W([@$rest, $e - 1], \@aa), R($mult))); }
    # - kappa_e
    if ($e == 0) { $res = vadd($res, vscale(W($rest, \@aa), R(-(2*$h - 2 + $n - 1)) * $mult)); }
    else { $res = vadd($res, vscale(W([@$rest, $e], \@aa), R(-$mult))); } }
  if ($an == 0) { for my $i (0..$#aa) { next if $aa[$i] < 1; my @c = @aa; $c[$i]--; $res = vadd($res, W(\@b, \@c)); } }
  $Wmemo{$key} = $res; return $res; }
1;
