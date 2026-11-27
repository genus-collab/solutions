# ref_loc5.pl -- referee check of the input (H) of Route D, Theorem D5 (no TR): virtual localization (Graber-Pandharipande) on the
# 5-fold Z = Tot(O(-1) + O(-1) + O + O -> P^1) = resolved conifold x C^2, with a torus acting ALSO on the base (alpha0 != 0):
#   tangent weights at p0: (alpha0; a1, a2, -e1, -e2), at p_oo: (-alpha0; a1 + alpha0, a2 + alpha0, -e1, -e2),
#   CY condition alpha0 + a1 + a2 = e1 + e2 (README of the project; Brini-Schuler, arXiv/main_DMJ.tex s. "The resolved conifold").
# Prediction (Brini-Schuler Prop. prop:ResConifold, proved by degeneration + Bryan-Pandharipande; hypothesis a_i = e_j for some i, j):
#   GW^Z_{g,d} = [degree 2g-2 part of] (1/d)(2 sinh(d e1/2))^-1 (2 sinh(d e2/2))^-1 = P_g d^(2g-3)/(e1 e2),
#   P_g = [u^(2g)] S(e1 u) S(e2 u),  S(y) = (y/2)/sinh(y/2),
# i.e. F_g - F_g(0) = P_g Li_{3-2g}(Q)/(e1 e2), whose only polar term at Q = 1 is (2g-3)! P_g/(e1 e2) t^(2-2g)  [this is (H)].
# Edge (degree d): prod_{k=1}^{d-1} (a1 + k alpha0/d)(a2 + k alpha0/d) / [ (-1)^d (d!)^2 alpha0^(2d)/d^(2d) ] / ((-e1)(-e2)), times 1/d.
# Stable vertex: e(T_p Z)^(val-1) Int prod_{w in T_p Z} Lambda^vee(w) prod_F 1/(t/d_F - psi_F); unstable (0,1): t/d; (0,2): e(T_p Z)/(t/d1 + t/d2).
# The Hodge integrals use the author's validated engine hodge_lib.pl (only for the vertex integrals; the localization code is the referee's).
# Usage: perl ref_loc5.pl   (runs the built-in list of tests)
use strict; use warnings; use FindBin;
{ my $d = "$FindBin::Bin/.."; require "$d/hodge_lib.pl"; }
require "$FindBin::Bin/ref_lib.pl";
sub S_coef { my $k = shift; # [y^(2k)] (y/2)/sinh(y/2) via 1/sinh series: sinh(y/2)/(y/2) = sum (y/2)^(2n)/(2n+1)!
  my @a = map { rq(1) / (rfact(2*$_ + 1) * rq(4)**$_) } 0..$k; my @b = (rq(1)); for my $n (1..$k) { my $s = rq(0); $s += $a[$_] * $b[$n-$_] for 1..$n; push @b, -$s; } return $b[$k]; }
sub P_g { my ($g, $e1, $e2) = @_; my $s = rq(0); for my $k (0..$g) { $s += S_coef($k) * S_coef($g - $k) * $e1**(2*$k) * $e2**(2*($g - $k)); } return $s; }
sub compositions { my ($n, $k) = @_; return ([]) if $k == 0 && $n == 0; return () if $k == 0 || $n < $k;
  my @out; for my $f (1..$n-$k+1) { push @out, [$f, @$_] for compositions($n - $f, $k - 1); } return @out; }
sub gen_assign { my ($tot, $k) = @_; return ([]) if $k == 0 && $tot == 0; return () if $k == 0;
  my @out; for my $f (0..$tot) { push @out, [$f, @$_] for gen_assign($tot - $f, $k - 1); } return @out; }
sub GW5 { my ($g, $d, $W) = @_; my ($al, $a1, $a2, $e1, $e2) = @$W;
  my @P = ([$al, $a1, $a2, -$e1, -$e2], [-$al, $a1 + $al, $a2 + $al, -$e1, -$e2]);
  my %VC;
  my $vert = sub { my ($p, $gv, @dF) = @_; my $key = join("|", $p, $gv, sort { $a <=> $b } @dF); return $VC{$key} if exists $VC{$key};
    my @w = @{$P[$p]}; my $t = $w[0]; my $eT = Q(1); $eT *= $_ for @w; my $val = scalar @dF; my $res;
    if ($gv == 0 && $val == 1) { $res = $t / $dF[0]; }
    elsif ($gv == 0 && $val == 2) { $res = $eT / ($t / $dF[0] + $t / $dF[1]); }
    else { my $cls = omega_class($gv, @w); my $N = 3*$gv - 3 + $val;
      my @ser = map { my $dd = $_; [ map { (Q($dd) / $t)**($_ + 1) } 0..$N ] } @dF;
      $res = $eT**($val - 1) * $eT**$gv * vertex($gv, $cls, @ser); }   # prod Lambda^vee(w) = (prod w)^g omega_class
    $VC{$key} = $res; return $res; };
  my $edge = sub { my $dd = shift; my $num = Q(1); $num *= ($a1 + Q($_) * $al / $dd) * ($a2 + Q($_) * $al / $dd) for 1..$dd-1;
    my $den = Q((-1)**$dd) * fact($dd)**2 * $al**(2*$dd) / Q($dd)**(2*$dd); return $num / $den / ($e1 * $e2) / $dd; };
  my $tot = Q(0);
  for my $E (1..$d) { for my $degs (compositions($d, $E)) { for my $n0 (1..$E) { for my $ni (1..$E) {
    my $h1 = $E - $n0 - $ni + 1; next if $h1 < 0 || $h1 > $g; my $np = $n0 * $ni;
    for my $code (0..$np**$E - 1) { my @pr; my $cc = $code; for (1..$E) { push @pr, $cc % $np; $cc = int($cc / $np); }
      my @v0 = (0) x $n0; my @vi = (0) x $ni; my @par = (0..$n0+$ni-1);
      my $find = sub { my $z = shift; $z = $par[$z] while $par[$z] != $z; return $z; };
      for my $e (0..$E-1) { my ($i, $j) = ($pr[$e] % $n0, int($pr[$e] / $n0)); $v0[$i]++; $vi[$j]++; my ($A, $B) = ($find->($i), $find->($n0 + $j)); $par[$A] = $B if $A != $B; }
      next if grep { $_ == 0 } (@v0, @vi); my $rt = $find->(0); next if grep { $find->($_) != $rt } 0..$n0+$ni-1;
      for my $gs (gen_assign($g - $h1, $n0 + $ni)) { my $cont = Q(1);
        $cont *= $edge->($degs->[$_]) for 0..$E-1;
        for my $i (0..$n0-1) { my @dF = map { $degs->[$_] } grep { $pr[$_] % $n0 == $i } 0..$E-1; $cont *= $vert->(0, $gs->[$i], @dF); last if $cont == 0; }
        next if $cont == 0;
        for my $j (0..$ni-1) { my @dF = map { $degs->[$_] } grep { int($pr[$_] / $n0) == $j } 0..$E-1; $cont *= $vert->(1, $gs->[$n0 + $j], @dF); last if $cont == 0; }
        $tot += $cont / (fact($n0) * fact($ni) * fact($E)); } } } } } }
  return $tot; }
my $bad = 0;
my ($e1, $e2) = (Q(2,3), Q(5,4));
my @tori = ( ["a1 = e1, alpha0 = 1",   [Q(1),   $e1, $e2 - 1]],
             ["a1 = e1, alpha0 = -3/7", [Q(-3,7), $e1, $e2 + Q(3,7)]],
             ["a2 = e2, alpha0 = 2",   [Q(2),  $e1 - 2, $e2]],
             ["generic T_CY (no a_i = e_j), alpha0 = 1, a1 = 1/3", [Q(1), Q(1,3), $e1 + $e2 - 1 - Q(1,3)]] );
my @cases = ([0,1], [0,2], [1,1], [1,2], [2,1], [2,2], [3,1]);
for my $T (@tori) { my ($name, $w) = @$T; my @W = ($w->[0], $w->[1], $w->[2], $e1, $e2);
  print "torus: $name (e1 = $e1, e2 = $e2); CY check alpha0 + a1 + a2 - e1 - e2 = ", $W[0] + $W[1] + $W[2] - $e1 - $e2, "\n";
  for my $cs (@cases) { my ($g, $d) = @$cs; next if $name =~ /generic/ && $g == 3;
    my $t0 = time; my $got = GW5($g, $d, \@W); my $pred = P_g($g, $e1, $e2) * Q($d)**(2*$g - 3) / ($e1 * $e2);
    my $okk = ($got == $pred); $bad++ unless $okk;
    printf "   (g,d) = (%d,%d): localization %s, prediction P_g d^(2g-3)/(e1 e2) = %s  %s [%d s]\n", $g, $d, $got, $pred, $okk ? "ok" : "MISMATCH", time - $t0; } }
print $bad ? "SOME MISMATCH ($bad)\n" : "ALL OK\n";
