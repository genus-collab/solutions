# slice_fitR.pl m1,m2,...,m_{r-1} [hmax] : the shortened-Norlund slice ansatz for r leaves R = (m_1, ..., m_{r-1}, j):
#   sigma^h Xi_h(m, j) = Sum_{s=0}^{M-1} prod_{l=s+1}^{N} (J - l) [lambda^{2h}] B^j C^{(m)}_s(lambda),
#   N = 2h-3+r, J = M + j, M = m_1 + ... + m_{r-1},  universal series C^{(m)}_s, [lambda^0] C_s = delta_{s0},
# determined genus by genus (h = 1..hmax); at genus h only s <= N (non-degenerate products) are used (ORD=2 convention).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use SerS;
sub R { Math::BigRat->new(shift) }
$| = 1;
my ($ms, $hmax) = @ARGV; $hmax //= 4; my @m = split /,/, $ms; my $r = scalar(@m) + 1; my $M = 0; $M += $_ for @m;
my @bp = (SerS::sone()); for my $p (1..$hmax) { push @bp, SerS::smul($bp[-1], SerS::beta()); }
sub fact { my $n = shift; my $x = R(1); $x *= $_ for 2..$n; $x }
sub BjAt { my ($k, $j) = @_; my $res = [R(0)]; for my $p (0..$k) { $res = SerS::padd($res, SerS::pscale(SerS::coeff($bp[$p], $k), R($j) ** $p / fact($p))); } return $res; }
sub nor { my ($s, $N, $J) = @_; my $x = R(1); $x *= ($J - $_) for ($s+1)..$N; return $x; }
my $S = $M >= 1 ? $M - 1 : 0;
my %C; $C{$_}{0} = [R($_ == 0 ? 1 : 0)] for 0..$S;
for my $h (1..$hmax) { Hodge::setup($h); %HodgePt::memo = (); my $N = 2*$h - 3 + $r;
  my @js = (0 .. 3*$h + $r + 2);
  my @res; for my $j (@js) { my $v = HodgePt::WP([], [@m, $j]); my @s = map { R(0) } 0..$h; for my $k (0..$#$v) { $s[$h - $k] += (($k % 2) ? -1 : 1) * $v->[$k]; }
    my $rr = [@s]; my $J = $M + $j;
    for my $s (0..$S) { for my $n (0..$h-1) { my $c = $C{$s}{$n} or next; my $t = SerS::pscale(SerS::pmul($c, BjAt($h - $n, $j)), nor($s, $N, $J)); $rr = SerS::padd($rr, SerS::pscale($t, R(-1))); } }
    push @res, $rr; }
  my $smax = $S < $N ? $S : $N; my $ok = 1; my @x;
  for my $sp (0..$h + 1) { my @rows = map { my $j = $js[$_]; [ (map { nor($_, $N, $M + $j) } 0..$smax), ($res[$_][$sp] // R(0)) ] } 0..$#js;
    my $nc = $smax + 1; my $row = 0; my @piv;
    for my $c (0..$nc-1) { my $p; for my $i ($row..$#rows) { if ($rows[$i][$c] != 0) { $p = $i; last } } next unless defined $p; @rows[$row, $p] = @rows[$p, $row];
      my $pv = $rows[$row][$c]; $_ = $_ / $pv for @{$rows[$row]}; for my $i (0..$#rows) { next if $i == $row || $rows[$i][$c] == 0; my $f = $rows[$i][$c]; $rows[$i][$_] -= $f * $rows[$row][$_] for 0..$nc; } push @piv, $c; $row++; }
    for my $i ($row..$#rows) { $ok = 0 if $rows[$i][$nc] != 0; }
    for my $k (0..$#piv) { $x[$piv[$k]][$sp] = $rows[$k][$nc]; } }
  for my $s (0..$S) { $C{$s}{$h} = [ map { ($x[$s] && defined $x[$s][$_]) ? $x[$s][$_] : R(0) } 0..$h+1 ]; }
  print "genus $h (R = (@m, j)): ", ($ok ? "consistent" : "ANSATZ FAILS"), "\n"; last unless $ok; }
for my $s (0..$S) { for my $n (sort { $a <=> $b } keys %{$C{$s}}) { my $p = $C{$s}{$n}; my @t = map { "($p->[$_])s^$_" } grep { ($p->[$_] // 0) != 0 } 0..$#$p; print "  C_$s [lambda^", 2*$n, "] = ", (@t ? join(' ', @t) : '0'), "\n"; } }
