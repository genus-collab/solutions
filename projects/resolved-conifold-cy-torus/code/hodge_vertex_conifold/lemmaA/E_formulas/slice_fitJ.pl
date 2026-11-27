# slice_fitJ.pl m [hmax] : test, on the slice R = (m, j), the "shortened Norlund" ansatz
#   sigma^h Xi_h(m, j) = Sum_{s=0}^{m-1} prod_{l=s+1}^{2h-1} (J - l) [lambda^{2h}] B^j C^{(m)}_s(lambda),   J = m + j,
# with universal series C^{(m)}_s (independent of h), [lambda^0] C_s = delta_{s0}; determined genus by genus (h = 1..hmax):
# at genus h the residual must lie in the span of the m polynomials prod_{l=s+1}^{2h-1}(J - l).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use SerS;
sub R { Math::BigRat->new(shift) }
$| = 1;
my ($m, $hmax) = @ARGV; $hmax //= 5;
my @bp = (SerS::sone()); for my $p (1..$hmax) { push @bp, SerS::smul($bp[-1], SerS::beta()); }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
sub BjAt { my ($k, $j) = @_; my $r = [R(0)]; for my $p (0..$k) { $r = SerS::padd($r, SerS::pscale(SerS::coeff($bp[$p], $k), R($j) ** $p / fact($p))); } return $r; }  # [lambda^{2k}] B^j
sub nor { my ($s, $h, $J) = @_; my $x = R(1); $x *= ($J - $_) for ($s+1)..(2*$h-1); return $x; }
my %C; $C{$_}{0} = [R($_ == 0 ? 1 : 0)] for 0..$m-1;
for my $h (1..$hmax) { Hodge::setup($h); %HodgePt::memo = ();
  my @js = (0 .. 3*$h + 3);
  my @res; for my $j (@js) { my $v = HodgePt::WP([], [$m, $j]); my @s = map { R(0) } 0..$h; for my $k (0..$#$v) { $s[$h - $k] += (($k % 2) ? -1 : 1) * $v->[$k]; }
    my $r = [@s]; my $J = $m + $j;
    for my $s (0..$m-1) { for my $n (0..$h-1) { my $c = $C{$s}{$n} or next; my $t = SerS::pscale(SerS::pmul($c, BjAt($h - $n, $j)), nor($s, $h, $J)); $r = SerS::padd($r, SerS::pscale($t, R(-1))); } }
    push @res, $r; }
  # solve residual(j) = Sum_s x_s nor(s,h,J) for each sigma power
  my $ok = 1; my @x;
  for my $sp (0..$h + 1) { my $ordc = $ENV{ORD} // 1; my $lim = $ordc == 1 ? $h : ($ordc == 2 ? 2*$h - 1 : 2*$h); my $smax = ($m - 1 < $lim ? $m - 1 : $lim);   # order condition: ORD=1: C_s = O(lambda^{2s}); ORD=2: s <= 2h-1 (no degenerate duplicates)
    my @rows = map { my $j = $js[$_]; [ (map { nor($_, $h, $m + $j) } 0..$smax), ($res[$_][$sp] // R(0)) ] } 0..$#js;
    my $nc = $smax + 1; my $row = 0; my @piv;
    for my $c (0..$nc-1) { my $p; for my $i ($row..$#rows) { if ($rows[$i][$c] != 0) { $p = $i; last } } next unless defined $p; @rows[$row, $p] = @rows[$p, $row];
      my $pv = $rows[$row][$c]; $_ = $_ / $pv for @{$rows[$row]}; for my $i (0..$#rows) { next if $i == $row || $rows[$i][$c] == 0; my $f = $rows[$i][$c]; $rows[$i][$_] -= $f * $rows[$row][$_] for 0..$nc; } push @piv, $c; $row++; }
    for my $i ($row..$#rows) { $ok = 0 if $rows[$i][$nc] != 0; }
    for my $k (0..$#piv) { $x[$piv[$k]][$sp] = $rows[$k][$nc]; } }
  for my $s (0..$m-1) { $C{$s}{$h} = [ map { ($x[$s] && defined $x[$s][$_]) ? $x[$s][$_] : R(0) } 0..$h+1 ]; }
  print "genus $h (slice m = $m): ", ($ok ? "consistent" : "ANSATZ FAILS"), "\n"; last unless $ok; }
for my $s (0..$m-1) { for my $n (sort { $a <=> $b } keys %{$C{$s}}) { my $p = $C{$s}{$n}; my @t = map { "($p->[$_])s^$_" } grep { ($p->[$_] // 0) != 0 } 0..$#$p; print "  C_$s [lambda^", 2*$n, "] = ", (@t ? join(' ', @t) : '0'), "\n"; } }
