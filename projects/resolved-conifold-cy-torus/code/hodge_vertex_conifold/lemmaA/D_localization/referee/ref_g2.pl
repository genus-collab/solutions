# ref_g2.pl -- referee's INDEPENDENT genus-2 evaluator (does not use hodge_lib.pl / xi_lib.pl):
#   <Omega prod psi^a>_2 (c-normalisation, mode-B class, x rational) is read from the Mathematica table lemmaA/data/psi_g2.txt (n <= 5);
#   larger n is reduced by the string and dilaton equations (Omega is pulled back from Mbar_2; in genus 2 every tuple with n >= 4
#   contains an entry 0 or 1, so the table suffices).
#   xi2(x, R)  = Xi^c_2(R) = Sum_m 1/m! < Omega prod 1/(c - j_i psi_i) prod_{k<=m} (-psi^2/(c - psi)) >_2
#   mv2(x, mu) = < Omega prod 1/(c - mu_i psi_i) >_2          (identical leaves are grouped)
# Uses: (a) cross-check against the q-tables (q_g2_r1, q_g2_r2 at s = -1); (b) the genus-2 bare-vertex test of Route D s.5 with the FULL
# leaf expansion from the closed form (ref_closedform.pl conventions) and with the truncated T_k form; (c) the Laguerre form of D2.3
# on cases the author did not use; (d) the k = 1 identity [t^1] F_h(rho q_t) = ((1-x)^2/2) Xi^c_h(2).
use strict; use warnings; use FindBin; require "$FindBin::Bin/ref_lib.pl";
my $data = "$FindBin::Bin/../../data";
my %TAB;   # "g|a1,a2,..." -> expression
for my $g (2) { open my $fh, '<', "$data/psi_g$g.txt" or die $!;
  while (my $l = <$fh>) { chomp $l; next unless $l =~ /\S/; my ($gg, $n, $as, $e) = split /\s*\|\s*/, $l, 4;
    $as =~ s/[{}\s]//g; my @a = length($as) ? split(/,/, $as) : (); $TAB{"$gg|" . join(",", sort { $a <=> $b } @a)} = $e; } close $fh; }
my %MEMO; our $CURX;
sub psi_int { my ($g, @a) = @_; @a = sort { $a <=> $b } @a; my $n = scalar @a; my $key = "$g|" . join(",", @a);
  my $mk = "$CURX|$key"; return $MEMO{$mk} if exists $MEMO{$mk};
  my $v; my $S = 0; $S += $_ for @a;
  if ($S > 3*$g - 3 + $n) { $v = rq(0); }
  elsif (exists $TAB{$key}) { $v = ev_expr($TAB{$key}, { x => $CURX }); }
  elsif ($n > 0 && $a[0] == 0) { my @r = @a[1..$n-1]; $v = rq(0); for my $j (0..$#r) { next if $r[$j] == 0; my @b = @r; $b[$j]--; $v += psi_int($g, @b); } }
  elsif ($n > 0 && $a[0] == 1) { my @r = @a[1..$n-1]; $v = (2*$g - 2 + scalar(@r)) * psi_int($g, @r); }
  else { die "psi_int: table entry missing for $key"; }
  $MEMO{$mk} = $v; return $v; }
# grouped enumeration: groups = list of [coef_sub(a) -> BigRat, multiplicity, min exponent]; sums Int prod coef(a_i) psi^a_i over exponent
# multisets per group (weight = multinomial / mult! = 1/prod(count!) for identical points, times the Aut-free labelled count)
sub grouped_integral { my ($g, $x, $groups) = @_; local $CURX = $x;
  my $n = 0; $n += $_->[1] for @$groups; my $dim = 3*$g - 3 + $n; my $tot = rq(0);
  my $rec; $rec = sub { my ($gi, $left, $exps, $w) = @_;
    if ($gi == @$groups) { $tot += $w * psi_int($g, @$exps); return; }
    my ($cf, $mult, $amin) = @{$groups->[$gi]};
    my $inner; $inner = sub { my ($k, $cap, $lft, $cur, $ww) = @_;
      if ($k == $mult) { my %cnt; $cnt{$_}++ for @$cur; my $lab = rfact($mult); $lab /= rfact($_) for values %cnt;
        $rec->($gi + 1, $lft, [@$exps, @$cur], $w * $ww * $lab); return; }
      my $top = $cap < $lft ? $cap : $lft;
      for my $a (reverse $amin..$top) { my $cc = $cf->($a); next if $cc == 0; $inner->($k + 1, $a, $lft - $a, [@$cur, $a], $ww * $cc); } };
    $inner->(0, $dim, $left, [], rq(1)); };
  $rec->(0, $dim, [], rq(1)); return $tot; }
sub leafgroups { my ($c, @R) = @_; my %m; $m{$_}++ for @R;
  return map { my $j = $_; [ sub { rq($j)**$_[0] / $c**($_[0] + 1) }, $m{$j}, 0 ] } sort { $a <=> $b } keys %m; }
sub xi2 { my ($x, @R) = @_; my $c = 1 + $x; my $r = scalar @R; my $tot = rq(0);
  for my $m (0..3*2 - 3 + $r) { my @G = leafgroups($c, @R);
    push @G, [ sub { -rq(1) / $c**($_[0] - 1) }, $m, 2 ] if $m > 0;
    $tot += grouped_integral(2, $x, \@G) / rfact($m); }     # labelled sum over the m dilaton points, divided by m!
  return $tot; }
sub mv2 { my ($x, @mu) = @_; my $c = 1 + $x; return grouped_integral(2, $x, [ leafgroups($c, @mu) ]); }
1 if caller;
return 1 if caller;

my $bad = 0; sub ok { my ($n, $a, $b) = @_; my $g = ($a == $b); $bad++ unless $g; printf "   %-70s %s\n", $n, $g ? "ok" : "MISMATCH ($a vs $b)"; }
# (a) against the q-tables
my %QT; for my $r (1, 2) { open my $fh, '<', "$data/q_g2_r$r.txt" or die $!; local $/; my $s = <$fh>; close $fh; my @p = split /\s*\|\s*/, $s; $QT{$r} = $p[-1]; }
for my $x (rq(2), rq(1,3), rq(-3,5), rq(5,2)) {
  print "x = $x\n";
  for my $R ([2], [3], [4], [5], [2,2], [2,3], [3,3], [2,4]) { my $r = scalar @$R; my %v = (x => $x, s => -1); $v{"j$_"} = $R->[$_-1] for 1..$r;
    ok("Xi^c_2(@$R): own evaluator vs q-table (s = -1)", xi2($x, @$R), ev_expr($QT{$r}, \%v)); }
  ok("Theorem C: Xi^c_2(1) = 0", xi2($x, 1), rq(0));
  ok("string identity Xi^c_2(1,3) = (3/c^2) Xi^c_2(3)", xi2($x, 1, 3), 3 / (1 + $x)**2 * xi2($x, 3));
}
# (b) genus-2 bare vertex: full leaf expansion (closed form) and truncated T_k form vs the direct values of leaf_structure.md, Step 2
require "$FindBin::Bin/ref_leafseries.pl";
for my $x (rq(2), rq(1,3), rq(-3,5), rq(5,2), rq(-7,3)) {
  my $c = 1 + $x; my $h = 2; print "x = $x: genus-2 bare vertex\n";
  my ($A, $C) = leaf_series($x, 4);
  my @Rs = ([2], [3], [2,2], [4], [2,3], [2,2,2]);
  my %X = map { join(",", @$_) => xi2($x, @$_) } @Rs;
  my %G = map { join(",", @$_) => G_R($A, $c, @$_) } @Rs;
  my $pre = (-$x * $c**3)**($h - 1);
  my @full = map { my $k = $_; my $s = rq(0); $s += $G{$_}[$k] * $X{$_} for keys %X; $pre * $s } (2, 3);
  my $T2 = $pre * (-$C->{3} * $X{"3"} + $C->{2}**2 / 2 * $X{"2,2"});
  my $T3 = $pre * ($C->{4} * $X{"4"} - $C->{2} * $C->{3} * $X{"2,3"} + $C->{2}**3 / 6 * $X{"2,2,2"});
  my $d0 = -(7 + 40*$x + 190*$x**2 - 280*$x**3 - 184*$x**4 - 280*$x**5 + 190*$x**6 + 40*$x**7 + 7*$x**8)/(69120*$x*$c**4);
  my $d1 = (1-$x)**2 * $x * (5 - 61*$x + 5*$x**2)/(3456*$c**4);
  ok("t^0: full leaf expansion = direct", $full[0], $d0); ok("t^0: truncated T_2 form = direct", $T2, $d0);
  ok("t^1: full leaf expansion = direct", $full[1], $d1); ok("t^1: truncated T_3 form = direct", $T3, $d1);
  ok("Xi_2(2) = 0 (so the D = k-2 residue at t^1 is harmless)", $X{"2"}, rq(0));
  # (d) the k = 1 identity
  ok("[t^1] G_(2) = (1-x)^2/2", $G{"2"}[1], (1 - $x)**2 / 2);
}
# (c) Laguerre form of D2.3 at genus 2, own evaluator, cases not used by the author
sub laguerre { my ($n, $al, $z) = @_; my $t = rq(0); for my $k (0..$n) { $t += rbinom($n + $al, $n - $k) * (-$z)**$k / rfact($k); } return $t; }
for my $x (rq(2), rq(-3,5)) { my $c = 1 + $x; my $tl = -$x / $c;
  for my $R ([4], [5], [2,3], [3,3], [2,2,2]) { my $h = 2; my $r = scalar @$R; my $J = 0; $J += $_ for @$R; my $a = 2*$h - 2 + $r; my $N0 = 3*$h - 3 + $r;
    # c = 1 normalisation: X^(c=1)(tau) = c^(3h-3+2r) X^c ; MV with l leaves: c^(3h-3+2l)
    my $sum = rq(0); for my $p (0..$N0) { my $phi = mv2($x, @$R, (1) x $p) * $c**(3*$h - 3 + 2*($r + $p)); $sum += rq((-1)**$p) * laguerre($N0 - $p, $a + $p, rq(-($J + $p))) * $phi / rfact($p); }
    my $xi = xi2($x, @$R) * $c**(3*$h - 3 + 2*$r);
    ok("Laguerre form, h = 2, R = (@$R), x = $x", $sum, $xi); } }
print $bad ? "SOME MISMATCH ($bad)\n" : "ALL OK\n";
