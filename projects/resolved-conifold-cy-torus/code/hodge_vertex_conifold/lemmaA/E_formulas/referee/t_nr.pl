# t_nr.pl hmax : independent test of the author's "Newton-regular" form NR (report 0(d)) for two leaves on the REFEREE data
# (ref_engine_h<h>_r2.txt, h = 1..hmax), written from scratch:
#   sigma^h Xi_h(j1,j2) = Sum_{(k,n')} (j1)_k (j2)_k prod_{l=n'+1}^{2h-1} (J-l) [L^h] ( S1 S2 B^J A_{k,n'} ),   L = lambda^2,
#   (k,n') in {(0,0)} u {k >= 1, k <= n' <= 2k-1},  A_{0,0} = 1,  A_{k,n'} = Sum_{i >= k} a_{k,n',i}(sigma) L^i  (h-independent).
# Genus by genus the new coefficients a_{k,n',h} (h(h+1)/2 unknown sigma-polynomials) are fitted to the residual on the full grid
# 0 <= j1 <= j2 <= 3h+1 (which determines a symmetric polynomial of degree <= 3h-1), separately for each power of sigma.
# Reports: number of equations, rank, consistency, and whether a_{k,n',h} = O(sigma^{2k-n'}) and the claimed values
# A_{1,1} = e^{-beta} - 1, leading A_{k,k} = (-sigma/8)^k/(k!(2k+1)) L^k come out.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
$| = 1; my $H = shift // 5;
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
# sigma-polynomials
sub padd { my ($a, $b) = @_; my $n = $#$a > $#$b ? $#$a : $#$b; return [ map { ($a->[$_] // Q(0)) + ($b->[$_] // Q(0)) } 0..$n ]; }
sub pmul { my ($a, $b) = @_; my @c = map { Q(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i] == 0; for my $j (0..$#$b) { $c[$i+$j] += $a->[$i] * $b->[$j] } } return \@c; }
sub psc { my ($a, $c) = @_; return [ map { $_ * $c } @$a ]; }
my $K = $H;
my @c = (Q(0)); for my $n (1..$K) { push @c, -RefHodge::bern(2*$n) / (2*$n * fact(2*$n)); }   # log S = Sum c_n x^{2n}
my @p = ([Q(2)], [Q(1)]); for my $n (2..2*$K) { push @p, padd($p[$n-1], pmul([Q(0), Q(1)], $p[$n-2])); }   # e1^n + e2^n, e1+e2 = 1, e1 e2 = -sigma
my @logSS = ([Q(0)], map { psc($p[2*$_], $c[$_]) } 1..$K);
my @beta = ([Q(0)], map { pmul([Q(0), Q(-1)], psc($p[2*$_-1], $c[$_])) } 1..$K);
sub sexp { my $a = shift; my @r = ([Q(1)]); for my $n (1..$K) { my $s = [Q(0)]; for my $m (1..$n) { $s = padd($s, psc(pmul($a->[$m], $r[$n-$m]), Q($m))); } $r[$n] = psc($s, Q(1)/$n); } return \@r; }
my %SB; sub SB { my $J = shift; return $SB{$J} //= sexp([ map { padd($logSS[$_], psc($beta[$_], Q($J))) } 0..$K ]); }   # S1 S2 B^J
# data sigma^h Xi_h(j1,j2) as sigma-polynomial
my %X;
for my $h (1..$H) { my @C; open my $fh, '<', "$FindBin::Bin/ref_engine_h${h}_r2.txt" or die "need ref_engine_h${h}_r2.txt"; while (my $line = <$fh>) { chomp $line; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4;
    for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; push @{$C[$k]}, [[split /,/, $m], Q($cf)]; } } close $fh;
  for my $j1 (0..3*$h + 1) { for my $j2 ($j1..3*$h + 1) { my $v = [map { Q(0) } 0..$h];
      for my $k (0..$h) { my $s = Q(0); for my $t (@{$C[$k] // []}) { my ($e, $cf) = @$t; $s += $cf * Q($j1) ** $e->[0] * Q($j2) ** $e->[1]; } $v->[$h - $k] += (($k % 2) ? -$s : $s); }
      $X{"$h|$j1|$j2"} = $v; } } }
my @lab = ([0, 0]); for my $k (1..$H) { for my $np ($k..2*$k-1) { push @lab, [$k, $np]; } }
my %A; $A{"0|0"}[0] = [Q(1)];   # A{"k|n'"}[i] = sigma-poly
sub basis { my ($k, $np, $h, $j1, $j2) = @_; my $J = $j1 + $j2; my $x = Q(1); for my $i (0..$k-1) { $x *= ($j1 - $i) * ($j2 - $i); } return $x if $x == 0; $x *= ($J - $_) for ($np+1)..(2*$h-1); return $x; }
for my $h (1..$H) { my @pts; for my $j1 (0..3*$h + 1) { for my $j2 ($j1..3*$h + 1) { push @pts, [$j1, $j2]; } }
  my @unk = grep { $_->[0] >= 1 && $_->[0] <= $h && $_->[1] <= 2*$h - 1 } @lab;
  # residual
  my @res; for my $pt (@pts) { my ($j1, $j2) = @$pt; my $r = [ @{$X{"$h|$j1|$j2"}} ]; my $sb = SB($j1 + $j2);
    for my $l (@lab) { my ($k, $np) = @$l; next if $np > 2*$h - 1; my $b = basis($k, $np, $h, $j1, $j2); next if $b == 0;
      for my $i (0..$h-1) { my $a = $A{"$k|$np"}[$i] or next; $r = padd($r, psc(pmul($a, $sb->[$h - $i]), -$b)); } }
    push @res, $r; }
  my $ok = 1; my $rank = 0; my $maxs = 0; $maxs = ($#$_ > $maxs ? $#$_ : $maxs) for @res;
  for my $s (0..$maxs) { my @rows = map { my $i = $_; [ (map { basis($_->[0], $_->[1], $h, @{$pts[$i]}) } @unk), ($res[$i][$s] // Q(0)) ] } 0..$#pts;
    my $nc = scalar @unk; my $row = 0; my @piv;
    for my $cc (0..$nc-1) { my $pv; for my $i ($row..$#rows) { if ($rows[$i][$cc] != 0) { $pv = $i; last } } next unless defined $pv; @rows[$row, $pv] = @rows[$pv, $row];
      my $d = $rows[$row][$cc]; $_ = $_ / $d for @{$rows[$row]}; for my $i (0..$#rows) { next if $i == $row || $rows[$i][$cc] == 0; my $f = $rows[$i][$cc]; $rows[$i][$_] -= $f * $rows[$row][$_] for 0..$nc; } push @piv, $cc; $row++; }
    for my $i ($row..$#rows) { $ok = 0 if $rows[$i][$nc] != 0; }
    $rank = $row if $row > $rank;
    for my $q (0..$#piv) { my ($k, $np) = @{$unk[$piv[$q]]}; my $x = $rows[$q][$nc]; next if $x == 0; my $cur = $A{"$k|$np"}[$h] // [Q(0)]; my @z = map { Q(0) } 0..$s; $z[$s] = $x; $A{"$k|$np"}[$h] = padd($cur, \@z); } }
  # sigma-order property
  my $ordbad = 0; for my $u (@unk) { my ($k, $np) = @$u; my $a = $A{"$k|$np"}[$h] or next; for my $s (0..$#$a) { $ordbad++ if $a->[$s] != 0 && $s < 2*$k - $np; } }
  printf "genus %d: %d grid points x %d sigma-powers, %d unknowns (rank %d): %s; coefficients violating O(sigma^{2k-n'}): %d\n", $h, scalar(@pts), $maxs + 1, scalar(@unk), $rank, ($ok ? "CONSISTENT" : "INCONSISTENT"), $ordbad;
  last unless $ok; }
# claimed closed values
my $emb = sexp([ map { psc($beta[$_], Q(-1)) } 0..$K ]);
my $ok11 = 1; for my $i (1..$H) { my $a = $A{"1|1"}[$i] // [Q(0)]; my $d = padd($a, psc($emb->[$i], Q(-1))); $ok11 = 0 if grep { $_ != 0 } @$d; }
print "A_{1,1} = e^{-beta} - 1 through L^$H: ", ($ok11 ? "yes" : "NO"), "\n";
for my $k (1..$H) { my $a = $A{"$k|$k"}[$k] // [Q(0)]; my $e = ((Q(-1) / 8) ** $k) / (fact($k) * (2*$k + 1)); my @ex = map { Q(0) } 0..$k; $ex[$k] = $e; my $d = padd($a, psc(\@ex, Q(-1)));
  printf "leading A_{%d,%d} = (-sigma/8)^%d/(%d! %d) L^%d: %s\n", $k, $k, $k, $k, 2*$k + 1, $k, ((grep { $_ != 0 } @$d) ? "NO" : "yes"); }
for my $l (@lab) { my ($k, $np) = @$l; next unless $k >= 1; for my $i ($k..$H) { my $a = $A{"$k|$np"}[$i] or next; next unless grep { $_ != 0 } @$a; printf "  A_{%d,%d} [L^%d] = %s\n", $k, $np, $i, join(' + ', map { "($a->[$_])s^$_" } grep { $a->[$_] != 0 } 0..$#$a); } }
