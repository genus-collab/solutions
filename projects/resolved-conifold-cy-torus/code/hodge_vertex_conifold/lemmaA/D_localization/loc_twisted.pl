# loc_twisted.pl -- localization on Mbar_{g,m}(X, d), X = Tot(O(-1)+O(-1) -> P^1), of the "dilaton-leaf" invariants
#   G^{S}_{g,d}(c) = sum_m 1/m! < prod_{k=1}^m gamma_S T_c(psi_k) >^X_{g,d},   T_c(psi) = -psi^2/(c - psi),
# with gamma_S = gamma_0 = [p0]/e(T_p0 X) (S = "0"), gamma_oo (S = "oo") or 1 = gamma_0 + gamma_oo (S = "both").
# Fixed-locus contributions (report.md, D2): a T-point lies on a stable vertex (K_c-dressing of the vertex) or at an
# unstable end of an edge (factor T_c(psi_x), psi_x = -t_p/d, replacing the Aut factor t_p/d).
# At c = w1 the p0-vertices are exactly the Lemma A vertices sum_m 1/m! Int Omega prod lam prod T (times weights).
# Usage: perl loc_twisted.pl g d S c w1 w2 w3
use strict; use warnings;
use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/xi_lib.pl"; }

sub E_edge { my ($d, $w1, $w2, $w3) = @_; my $num = Q(1);
  for my $k (1..$d-1) { $num *= ($w2 + Q($k) * $w1 / $d) * ($w3 + Q($k) * $w1 / $d); }
  my $den = Q((-1)**$d) * fact($d)**2 * $w1**(2*$d) / Q($d)**(2*$d);
  return $num / $den; }
sub Tc { my ($c, $psi) = @_; return -$psi**2 / ($c - $psi); }

my %VC;
sub V_vertex {   # $p = 0 or 1 (p_oo); $tw = 1 if T-points allowed at this fixed point
  my ($p, $g, $W, $c, $tw, @dF) = @_; my ($w1, $w2, $w3) = @$W;
  my ($t, $f2, $f3) = $p == 0 ? ($w1, $w2, $w3) : (-$w1, $w2 + $w1, $w3 + $w1);
  my $val = scalar @dF;
  my $key = join("|", $p, $g, $tw, sort { $a <=> $b } @dF);
  return $VC{$key} if exists $VC{$key};
  my $res;
  if ($g == 0 && $val == 1) { $res = $t / $dF[0]; $res += Tc($c, -$t / $dF[0]) if $tw; }
  elsif ($g == 0 && $val == 2) { $res = $t * $f2 * $f3 / ($t / $dF[0] + $t / $dF[1]); }
  else {
    my $cls = omega_class($g, $t, $f2, $f3);
    my $N = 3*$g - 3 + $val + 2;
    # leaves 1/(t/d - psi) = (d/t) * 1/(1 - (d/t) psi): use the component machinery with scales d/t
    my @R = map { Q($_) / $t } @dF; my $pref = Q(1); $pref *= $_ for @R;
    my $I;
    if ($tw) { my $Nt = 3*$g - 3 + $val + 12; my $T = ser_T($c, $Nt); $I = components($g, $T, "Tc=$c", @R); }
    else { $I = components($g, undef, "none", @R); }
    $res = ($t * $f2 * $f3)**($val - 1 + $g) * $pref * contract($cls, $I);
  }
  $VC{$key} = $res; return $res;
}
sub compositions { my ($n, $k) = @_; return ([]) if $k == 0 && $n == 0; return () if $k == 0 || $n < $k;
  my @out; for my $first (1..$n-$k+1) { push @out, [$first, @$_] for compositions($n - $first, $k - 1); } return @out; }
sub genus_assign { my ($tot, $k) = @_; return ([]) if $k == 0 && $tot == 0; return () if $k == 0;
  my @out; for my $first (0..$tot) { push @out, [$first, @$_] for genus_assign($tot - $first, $k - 1); } return @out; }
sub G_gd {
  my ($g, $d, $S, $c, $W) = @_; my $total = Q(0);
  my ($tw0, $twI) = ($S eq "0") ? (1, 0) : ($S eq "oo") ? (0, 1) : ($S eq "both") ? (1, 1) : (0, 0);
  # degree 0: a single vertex with no edges (only if g >= 2): (t f2 f3)^{g-1} Int Omega_t K_c
  if ($d == 0) { for my $p (0, 1) { my $tw = $p == 0 ? $tw0 : $twI; next unless $tw || 1;
      my ($w1, $w2, $w3) = @$W; my ($t, $f2, $f3) = $p == 0 ? ($w1, $w2, $w3) : (-$w1, $w2 + $w1, $w3 + $w1);
      my $cls = omega_class($g, $t, $f2, $f3); my $I = $tw ? components($g, ser_T($c, 3*$g + 12), "Tc=$c") : components($g, undef, "none");
      $total += ($t * $f2 * $f3)**($g - 1) * contract($cls, $I); } return $total; }
  for my $E (1..$d) { for my $degs (compositions($d, $E)) {
    for my $n0 (1..$E) { for my $ninf (1..$E) {
      next if $n0 + $ninf > $E + 1; my $h1 = $E - $n0 - $ninf + 1; next if $h1 < 0 || $h1 > $g;
      my $npairs = $n0 * $ninf;
      for my $code (0..$npairs**$E - 1) {
        my @pairs; my $cc = $code; for my $e (0..$E-1) { push @pairs, $cc % $npairs; $cc = int($cc / $npairs); }
        my @val0 = (0) x $n0; my @valI = (0) x $ninf; my @par = (0..$n0+$ninf-1);
        my $find = sub { my $x = shift; while ($par[$x] != $x) { $x = $par[$x]; } return $x; };
        for my $e (0..$E-1) { my $i = $pairs[$e] % $n0; my $j = int($pairs[$e] / $n0); $val0[$i]++; $valI[$j]++;
          my ($a, $b) = ($find->($i), $find->($n0 + $j)); $par[$a] = $b if $a != $b; }
        next if grep { $_ == 0 } (@val0, @valI);
        my $root = $find->(0); my $conn = 1; for my $x (0..$n0+$ninf-1) { $conn = 0 if $find->($x) != $root; } next unless $conn;
        for my $gs (genus_assign($g - $h1, $n0 + $ninf)) {
          my $cont = Q(1);
          for my $e (0..$E-1) { $cont *= E_edge($degs->[$e], @$W) / $degs->[$e]; }
          for my $i (0..$n0-1) { my @dF; for my $e (0..$E-1) { push @dF, $degs->[$e] if $pairs[$e] % $n0 == $i; }
            $cont *= V_vertex(0, $gs->[$i], $W, $c, $tw0, @dF); last if $cont == 0; }
          next if $cont == 0;
          for my $j (0..$ninf-1) { my @dF; for my $e (0..$E-1) { push @dF, $degs->[$e] if int($pairs[$e] / $n0) == $j; }
            $cont *= V_vertex(1, $gs->[$n0 + $j], $W, $c, $twI, @dF); last if $cont == 0; }
          $total += $cont / (fact($n0) * fact($ninf) * fact($E));
        } } } } } }
  return $total;
}
if (!caller) {
  my ($g, $d, $S, $c, @w) = @ARGV; my @W = map { Q($_) } @w;
  print "G^{$S}_{$g,$d}(c = $c; w = @w) = ", G_gd($g, $d, $S, Q($c), \@W), "\n";
}
1;
