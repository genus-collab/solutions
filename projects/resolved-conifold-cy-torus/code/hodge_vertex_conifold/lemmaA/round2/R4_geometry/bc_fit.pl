# bc_fit.pl gmax cmax : y-Taylor coefficients B^(c)_g(j) = [y^c] Xi^(0)_g(j,y) of the two-leaf Hurwitz polynomial, from Route B's
# two-leaf REC expanded in y (no Hodge engine):
#   (2h+1)[(j-2h) B^(c)_h + B^(c-1)_h] = Sum_{h''<h} c_{h+1-h''} Sum_{a=0}^{c} (1/a!) d^a/dx^a (x)_(2h+2-2h'') |_{x=j-2h''} B^(c-a)_{h''}
#                                       + Sum_{m1+m2=h} T_{m1}(j) [y^c] T_{m2}(y),
# B^(0)_g = T_g(j)/j, genus 0: B^(c)_0 = (-1)^c / j^{c+1}.  Then tests the ansatz (for integers 1 <= j <= 2g-1)
#   B^(c)_g(j) = Sum_{g1=0}^{g-1} A^(c)_{g-g1} T_{g1}(j)/j   with universal constants A^(c)_m (fitted at genus m, tested at higher genus).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use PJ;
sub Q { Math::BigRat->new(shift) }
$| = 1; my ($gmax, $cmax) = @ARGV; $gmax //= 7; $cmax //= 3; $PJ::N = 2*$gmax + 4;
my $f = PJ::logS(); my $tf1 = PJ::stder($f); my $logH = PJ::sadd($f, $tf1, Q(-1));
my $HJ = PJ::sexp(PJ::sscalp($logH, PJ::pj()));
my %T; for my $n (0..$gmax) { $T{$n} = PJ::pmul(PJ::ffj(-1, 2*$n), $HJ->[2*$n]); }
sub c { my $k = shift; return Q(2) / PJ::fact(2*$k) }
# a-th derivative of (x)_(n) at x = j - s, as polynomial in j: d^a/dx^a applied to ffj(-s, n) in j (same thing)
sub dpow { my ($p, $a) = @_; my $q = $p; $q = PJ::pder($q) for 1..$a; return $q }
# rational functions in j: represent B^(c)_g as polynomial for g >= 1 (they are), genus 0 handled separately by evaluation
my %B;   # B{g}{c} polynomial in j (g >= 1)
for my $h (1..$gmax) {
  for my $c (0..$cmax) {
    if ($c == 0) { my ($q, $r) = PJ::pdiv_lin($T{$h}, 0); die if $r != 0; $B{$h}{0} = $q; next; }
    # RHS polynomial; genus-0 one-block term: c_{h+1} (j+y)_(2h+2)/(j+y) = c_{h+1} (j+y-1)_(2h+1) -> [y^c] = (1/c!) d^c/dj^c (j-1)_(2h+1)
    my $rhs = PJ::pscale(dpow(PJ::ffj(-1, 2*$h + 1), $c), c($h + 1) / PJ::fact($c));
    for my $hp (1..$h-1) { my $k = $h + 1 - $hp; my $ffp = PJ::ffj(-2*$hp, 2*$k);
      for my $a (0..$c) { $rhs = PJ::padd($rhs, PJ::pscale(PJ::pmul(dpow($ffp, $a), $B{$hp}{$c - $a}), c($k) / PJ::fact($a))); } }
    # two-block: Sum_{m1+m2=h} T_{m1}(j) [y^c] T_{m2}(y)
    for my $m2 (0..$h) { my $cy = $T{$m2}->[$c] // Q(0); next if $cy == 0; $rhs = PJ::padd($rhs, PJ::pscale($T{$h - $m2}, $cy)); }
    $rhs = PJ::psub($rhs, PJ::pscale($B{$h}{$c - 1}, Q(2*$h + 1)));
    my ($q, $r) = PJ::pdiv_lin($rhs, 2*$h); die "genus $h c $c: not divisible (rem $r)" if $r != 0;
    $B{$h}{$c} = PJ::pscale($q, Q(1) / (2*$h + 1)); } }
# validation against the two-leaf engine tables
for my $h (1..($gmax < 7 ? $gmax : 7)) { my $fn = "$FindBin::Bin/nr_test/engine_h${h}_r2.txt"; next unless -e $fn; open my $fh, '<', $fn or die; my %tab;
  while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; next unless $k == 0; for my $t (split /\s*;\s*/, $rest) { my ($m, $v) = split /:/, $t; my ($a, $b) = split /,/, $m; $tab{$b}[$a] = Q($v); } } close $fh;
  my $bad = 0; for my $c (0..$cmax) { my $p = $B{$h}{$c}; my $t = $tab{$c} // []; my $n = $#$p > $#$t ? $#$p : $#$t; for my $a (0..$n) { $bad++ if ($p->[$a] // Q(0)) != ($t->[$a] // Q(0)); } }
  print "genus $h: REC y-Taylor coefficients c <= $cmax vs two-leaf table: ", ($bad ? "$bad MISMATCHES" : "equal"), "\n"; }
# ansatz fit: at genus g, unknown A^(c)_g appears with T_0(j)/j = 1/j; lower A's known from lower genera
for my $c (1..$cmax) { my %A; my $ok = 1; my @rep;
  for my $g (1..$gmax) {
    # determine A^(c)_g from j = 1:  T_{g1}(1) = 0 for g1 >= 1, T_0 = 1  =>  B^(c)_g(1) = A^(c)_g
    $A{$g} = PJ::peval($B{$g}{$c}, Q(1));
    my @good; for my $j (1..2*$g) { my $jq = Q($j); my $s = Q(0); for my $g1 (0..$g-1) { $s += $A{$g - $g1} * PJ::peval($T{$g1}, $jq) / $jq; }
      push @good, $j if $s == PJ::peval($B{$g}{$c}, $jq); }
    push @rep, "g=$g: j=" . join(',', @good); }
  print "c=$c: ansatz B^(c)_g = Sum A^(c)_{g-g1} T_{g1}/j holds at  ", join(' | ', @rep), "\n";
  print "      A^(c)_m (m=1..$gmax): ", join(', ', map { $A{$_} } 1..$gmax), "\n"; }
