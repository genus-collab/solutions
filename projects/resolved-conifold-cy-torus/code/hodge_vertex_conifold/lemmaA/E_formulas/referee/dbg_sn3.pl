# t_sn.pl hmax m1[,m2,...] : independent test of the author's one-leaf universality (SN, report 0(c)) on the REFEREE data:
#   sigma^h Xi_h(m, j) = Sum_{s=0}^{max(M-1,0)} prod_{l=s+1}^{N} (M + j - l) [L^h] ( B^j C_s(L) ),   M = Sum m_i, N = 2h-3+r, L = lambda^2,
# C_s universal (h-independent), [L^0] C_s = delta_{s0}, C_s = O(L^s) (the author's ORD = 1).  Genus by genus the new coefficients
# [L^h] C_s (s <= min(M-1, h)) are fitted per power of sigma on the grid j = 0..3h+3 (determines the slice polynomial), and the
# remaining h conditions per sigma-power are the test.  Data: ref_engine_h<h>_r<r>.txt (or the compact ref_mset_h<h>_r<r>.txt).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($H, $ms) = @ARGV; my @m = split /,/, $ms; my $r = scalar(@m) + 1; my $M = 0; $M += $_ for @m;
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
sub padd { my ($a, $b) = @_; my $n = $#$a > $#$b ? $#$a : $#$b; return [ map { ($a->[$_] // Q(0)) + ($b->[$_] // Q(0)) } 0..$n ]; }
sub pmul { my ($a, $b) = @_; my @c = map { Q(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { next if $a->[$i] == 0; for my $j (0..$#$b) { $c[$i+$j] += $a->[$i] * $b->[$j] } } return \@c; }
sub psc { my ($a, $c) = @_; return [ map { $_ * $c } @$a ]; }
my $K = $H;
my @c = (Q(0)); for my $n (1..$K) { push @c, -RefHodge::bern(2*$n) / (2*$n * fact(2*$n)); }
my @p = ([Q(2)], [Q(1)]); for my $n (2..2*$K) { push @p, padd($p[$n-1], pmul([Q(0), Q(1)], $p[$n-2])); }
my @beta = ([Q(0)], map { pmul([Q(0), Q(-1)], psc($p[2*$_-1], $c[$_])) } 1..$K);
sub sexp { my $a = shift; my @r = ([Q(1)]); for my $n (1..$K) { my $s = [Q(0)]; for my $mm (1..$n) { $s = padd($s, psc(pmul($a->[$mm], $r[$n-$mm]), Q($mm))); } $r[$n] = psc($s, Q(1)/$n); } return \@r; }
my %BJ; sub BJ { my $j = shift; return $BJ{$j} //= sexp([ map { psc($beta[$_], Q($j)) } 0..$K ]); }
sub value_at { my ($h, $pt) = @_;   # sigma^h Xi_h(pt) as sigma-polynomial, from the referee files
  our %CACHE; my $C = $CACHE{$h} //= do { my @C; my $f = "$FindBin::Bin/ref_engine_h${h}_r${r}.txt"; my $compact = 0;
      unless (-e $f) { $f = "$FindBin::Bin/ref_mset_h${h}_r${r}.txt"; $compact = 1; }
      open my $fh, '<', $f or die "no data for h=$h r=$r"; while (my $line = <$fh>) { chomp $line; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4;
        for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($mm, $cf) = split /:/, $t; push @{$C[$k]}, [[split /,/, $mm], Q($cf)]; } } close $fh; [\@C, $compact] };
  my ($CC, $compact) = @$C; my $v = [map { Q(0) } 0..$h];
  for my $k (0..$h) { my $s = Q(0); for my $t (@{$CC->[$k] // []}) { my ($e, $cf) = @$t;
      if ($compact) { $s += $cf * RefHodge::mA([map { Q($_) } @$pt], $e); } else { my $x = $cf->copy; for my $i (0..$r-1) { $x *= Q($pt->[$i]) ** $e->[$i] if $e->[$i]; } $s += $x; } }
    $v->[$h - $k] += (($k % 2) ? -$s : $s); }
  return $v; }
my $smaxM = $M - 1 > 0 ? $M - 1 : 0;
my %C; $C{0}[0] = [Q(1)];
for my $h (1..$H) { my $N = 2*$h - 3 + $r; my @js = (0..3*$h + 3);
  my @unk = grep { $_ <= $h } 0..$smaxM;
  my @res; for my $j (@js) { my $v = value_at($h, [@m, $j]); my $bj = BJ($j);
    for my $s (0..$smaxM) { my $nor = Q(1); $nor *= ($M + $j - $_) for ($s+1)..$N; next if $nor == 0;
      for my $i (0..$h-1) { my $a = $C{$s}[$i] or next; $v = padd($v, psc(pmul($a, $bj->[$h - $i]), -$nor)); } }
    push @res, $v; }
  if ($ENV{DBG}) { for my $q (0..$#js) { print "h=$h j=$js[$q] res=", join(" ", @{$res[$q]}), "
"; } }
  my $ok = 1; my $maxs = 0; $maxs = ($#$_ > $maxs ? $#$_ : $maxs) for @res; my $neq = 0;
  for my $sp (0..$maxs) { my @rows = map { my $jj = $js[$_]; [ (map { my $s = $_; my $nor = Q(1); $nor *= ($M + $jj - $_) for ($s+1)..$N; $nor } @unk), ($res[$_][$sp] // Q(0)) ] } 0..$#js;
    my $nc = scalar @unk; my $row = 0; my @piv;
    for my $cc (0..$nc-1) { my $pv; for my $i ($row..$#rows) { if ($rows[$i][$cc] != 0) { $pv = $i; last } } next unless defined $pv; @rows[$row, $pv] = @rows[$pv, $row];
      my $d = $rows[$row][$cc]; $_ = $_ / $d for @{$rows[$row]}; for my $i (0..$#rows) { next if $i == $row || $rows[$i][$cc] == 0; my $f = $rows[$i][$cc]; $rows[$i][$_] -= $f * $rows[$row][$_] for 0..$nc; } push @piv, $cc; $row++; }
    for my $i ($row..$#rows) { $ok = 0 if $rows[$i][$nc] != 0; }
    for my $q (0..$#piv) { my $s = $unk[$piv[$q]]; my $x = $rows[$q][$nc]; next if $x == 0; my @z = map { Q(0) } 0..$sp; $z[$sp] = $x; $C{$s}[$h] = padd($C{$s}[$h] // [Q(0)], \@z); } }
  printf "SN slice m = (%s), genus %d: %d grid values, %d unknowns per sigma-power: %s\n", $ms, $h, scalar(@js), scalar(@unk), ($ok ? "consistent" : "ANSATZ FAILS"); last unless $ok; }
for my $s (0..$smaxM) { for my $i (0..$H) { my $a = $C{$s}[$i] or next; next unless grep { $_ != 0 } @$a; printf "  C_%d [L^%d] = %s\n", $s, $i, join(' ', map { "($a->[$_])s^$_" } grep { $a->[$_] != 0 } 0..$#$a); } }
