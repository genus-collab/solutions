#!/usr/bin/perl
# The leading-order mKP relation M_{G,eps}(R) for one-part targets R = (j), with the genus >= 2 values taken from the
# Mathematica tables lemmaA/data/q_g<g>_r<r>.txt (Mumford's formula; independent of Route E's engine):
#   Xi_g(j1..jr) (c = 1) = c^(3g-3+2r) Q_g(s = -1; j; x),  x = -tau/(1+tau), c = 1 + x  (brief s.1).
# Genus 0: closed forms; genus 1: Route E tables engine_h1_r*.txt (exact).  First prints a cross-check of the
# q-table values against the engine values (XiData) at a few points.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use mvlib (); use XiData;
$| = 1;
sub R { Math::BigRat->new(@_) }
my $DDIR = "$FindBin::Bin/../../data";
my %Q;   # "g|r" -> AST
for my $g (2..4) { for my $r (1..2) { my $f = "$DDIR/q_g${g}_r${r}.txt"; open(my $fh, '<', $f) or die "$f"; local $/; my $t = <$fh>; close $fh;
  $t =~ s/\r//g; my @p = split /\s*\|\s*/, $t, 3; $Q{"$g|$r"} = mvlib::parse_expr($p[2]); } }
my %qc;
sub Xq { my ($g, $T, $tau) = @_; my @T = sort { $a <=> $b } @$T; my $r = scalar @T; my $key = "$g|@T|$tau"; return $qc{$key} if exists $qc{$key};
  my $x = -$tau / (1 + $tau); my $c = 1 + $x; my %num = (x => $x, s => R(-1)); $num{"j$_"} = R($T[$_-1]) for 1..$r;
  my $p = mvlib::eval_poly($Q{"$g|$r"}, [], \%num); my $v = mvlib::pconstval($p) * ($c ** (3*$g - 3 + 2*$r)); $qc{$key} = $v; return $v; }
sub X { my ($g, $T, $tau) = @_; return XiData::xi0($T) if $g == 0; return XiData::xi_at(1, $T, $tau) if $g == 1;
  die "no q data for genus $g r=" . scalar(@$T) unless $g <= 4 && @$T <= 2; return Xq($g, $T, $tau); }
# cross-check of the data sources
my $nx = 0; my $bx = 0;
for my $tau (R("2/7"), R("-5/3")) { for my $g (2..4) { for my $T ([2], [3], [5], [8], [2,2], [2,5], [3,4], [2,9]) {
  my $a = Xq($g, $T, $tau); my $b = XiData::xi_at($g, $T, $tau); $nx++; $bx++ if $a != $b; } } }
print "data cross-check q-tables vs engine: $nx values, $bx mismatches\n";
sub fact { my $n = shift; my $f = R(1); $f *= $_ for 2..$n; $f }
sub ff { my ($x, $k) = @_; my $p = R(1); for my $i (0..$k-1) { $p *= ($x - $i); } $p }
sub Pi { my ($k, $tau) = @_; my $p = R(1); for my $a (1..$k-1) { $p *= ($k*$tau + $a); } $p / fact($k-1) }
sub Phi { my ($g, $k, $T, $tau, $sig) = @_;
  if (!@$T) { return R(0) if $g > 0; return $tau if $k == 1; return $tau*$tau/2 if $k == 2; return R(0); }
  my $J = 0; $J += $_ for @$T; my $D = $J - scalar(@$T); my $f = ff(R($D - 2*$g + 2), $k); return R(0) if $f == 0;
  return ($sig ** $g) * ($tau ** $k) * $J * $f * X($g, $T, $tau) / fact($k); }
sub rel { my ($G, $j, $tau, $par) = @_; my $sig = $tau*($tau+1); my $P2 = Pi(2,$tau); my $n = 2*$G + $par; my $J = $j; my $D = $j - 1;
  my $tot = R(0);
  for my $g (0..$G+1) { my $k = $n + 1 - 2*$g; next if $k < 1; my $f = ff(R($j + 2 - 2 - 2*$g + 2), $k); next if $f == 0;
    $tot += -2*$P2*($sig**$g)*($tau**$k)*$f*X($g, [2,$j], $tau)/fact($k); }
  $tot += 2*($sig**($G+1))*$J*($J+1)*X($G, [$j], $tau) if $par == 0;
  for my $g (0..$G) { my $k = $n - 2*$g; next if $k < 1; my $f = ff(R($D - 2*$g + 2), $k); next if $f == 0;
    $tot += ($sig**($g+1))*($tau**$k)*$J*($J+1)*$f*X($g, [$j], $tau)/fact($k); }
  for my $s ([[ $j ], []], [[], [ $j ]]) { my ($T1, $T2) = @$s;
    for my $g1 (0..int(($n+1)/2)) { for my $k1 (1..$n+1-2*$g1) { my $rest = $n + 2 - 2*$g1 - $k1; my $f1 = Phi($g1,$k1,$T1,$tau,$sig); next if $f1 == 0;
        for my $g2 (0..int(($rest-1)/2)) { my $k2 = $rest - 2*$g2; next if $k2 < 1; my $f2 = Phi($g2,$k2,$T2,$tau,$sig); next if $f2 == 0; $tot += $f1*$f2; } } } }
  return $tot; }
my ($n, $bad) = (0, 0);
for my $G (2..4) { for my $j (2..11) { for my $par (0, 1) { my @res = map { rel($G, $j, R($_), $par) } ("2/7", "3/2", "-5/3");
  my $ok = !(grep { $_ != 0 } @res); $n++; $bad++ unless $ok;
  printf "G=%d R=(%d) %s [R+2: D=%d excess=%d] %s\n", $G, $j, ($par ? "ODD " : "EVEN"), $j, 2*$G-2-$j, ($ok ? "OK" : "FAIL @res"); } } }
printf "SUMMARY (q-table data): %d relations, %d failures\n", $n, $bad;
