# loc_conifold.pl -- virtual localization for X = Tot(O(-1)+O(-1) -> P^1) with a torus T acting with weights
#   p0: (w1, w2, w3) (w1 = tangent weight of P^1, w2, w3 fibre weights),  p_oo: (-w1, w2 + w1, w3 + w1).
# Graber-Pandharipande formula (conventions in report.md, section D2):
#   N_{g,d} = sum_Gamma 1/|Aut Gamma| prod_e E(d_e)/d_e prod_v V(v),
#   E(d) = prod_{k=1}^{d-1} (w2 + k w1/d)(w3 + k w1/d) * d^{2d} / ((-1)^d (d!)^2 w1^{2d}),
#   V(v) = (t f2 f3)^{val-1} Int_{Mbar_{g,val}} Lambda^vee(t) Lambda^vee(f2) Lambda^vee(f3) prod_F 1/(t/d_F - psi_F),
#   unstable: (0,1): t/d ; (0,2): (t f2 f3)/(t/d1 + t/d2),
# where (t, f2, f3) are the weights at the vertex's fixed point.  Lambda^vee(u) = sum (-1)^i lambda_i u^{g-i}.
# Usage: perl loc_conifold.pl g d w1 w2 w3   (rationals; for the CY torus take w3 = -w1-w2)
use strict; use warnings;
use File::Basename qw(dirname); { my $d = dirname(__FILE__); require "$d/hodge_lib.pl"; }

sub E_edge { my ($d, $w1, $w2, $w3) = @_; my $num = Q(1);
  for my $k (1..$d-1) { $num *= ($w2 + Q($k) * $w1 / $d) * ($w3 + Q($k) * $w1 / $d); }
  my $den = Q((-1)**$d) * fact($d)**2 * $w1**(2*$d) / Q($d)**(2*$d);
  return $num / $den; }

my %VC;
sub V_vertex {   # $p = 0 or 1 (p_oo), genus, list of flag degrees
  my ($p, $g, $W, @dF) = @_; my ($w1, $w2, $w3) = @$W;
  my ($t, $f2, $f3) = $p == 0 ? ($w1, $w2, $w3) : (-$w1, $w2 + $w1, $w3 + $w1);
  my $val = scalar @dF;
  my $key = join("|", $p, $g, sort { $a <=> $b } @dF);
  return $VC{$key} if exists $VC{$key};
  my $res;
  if ($g == 0 && $val == 1) { $res = $t / $dF[0]; }
  elsif ($g == 0 && $val == 2) { $res = $t * $f2 * $f3 / ($t / $dF[0] + $t / $dF[1]); }
  else {
    # Lambda^vee(t) Lambda^vee(f2) Lambda^vee(f3) = (t f2 f3)^g prod c_{1/w}(E^vee)
    my $cls = omega_class($g, $t, $f2, $f3);
    my $N = 3*$g - 3 + $val;
    my @ser; for my $dd (@dF) { # 1/(t/d - psi) = (d/t) sum (d psi/t)^a
      my @s; my $c = Q($dd) / $t; for my $a (0..$N) { push @s, $c * (Q($dd) / $t)**$a; } push @ser, \@s; }
    $res = ($t * $f2 * $f3)**($val - 1) * ($t * $f2 * $f3)**$g * vertex($g, $cls, @ser);
  }
  $VC{$key} = $res; return $res;
}

# enumerate labelled bipartite graphs; returns sum over iso classes of Cont/|Aut|
sub compositions { my ($n, $k) = @_; return ([]) if $k == 0 && $n == 0; return () if $k == 0 || $n < $k;
  my @out; for my $first (1..$n-$k+1) { push @out, [$first, @$_] for compositions($n - $first, $k - 1); } return @out; }
sub genus_assign { my ($tot, $k) = @_; return ([]) if $k == 0 && $tot == 0; return () if $k == 0;
  my @out; for my $first (0..$tot) { push @out, [$first, @$_] for genus_assign($tot - $first, $k - 1); } return @out; }

sub N_gd {
  my ($g, $d, $W) = @_; my $total = Q(0);
  for my $E (1..$d) {
    for my $degs (compositions($d, $E)) {
      for my $n0 (1..$E) { for my $ninf (1..$E) {
        next if $n0 + $ninf > $E + 1;
        my $h1 = $E - $n0 - $ninf + 1; next if $h1 < 0 || $h1 > $g;
        my $npairs = $n0 * $ninf;
        # assignments of each edge to a pair (i, j)
        my $nassign = $npairs**$E;
        for my $code (0..$nassign-1) {
          my @pairs; my $c = $code; for my $e (0..$E-1) { push @pairs, $c % $npairs; $c = int($c / $npairs); }
          # connectivity and valence >= 1
          my @val0 = (0) x $n0; my @valI = (0) x $ninf;
          my @par = (0..$n0+$ninf-1);
          my $find; $find = sub { my $x = shift; while ($par[$x] != $x) { $x = $par[$x]; } return $x; };
          for my $e (0..$E-1) { my $i = $pairs[$e] % $n0; my $j = int($pairs[$e] / $n0); $val0[$i]++; $valI[$j]++;
            my ($a, $b) = ($find->($i), $find->($n0 + $j)); $par[$a] = $b if $a != $b; }
          next if grep { $_ == 0 } (@val0, @valI);
          my $root = $find->(0); my $conn = 1; for my $x (0..$n0+$ninf-1) { $conn = 0 if $find->($x) != $root; } next unless $conn;
          for my $gs (genus_assign($g - $h1, $n0 + $ninf)) {
            my $cont = Q(1);
            for my $e (0..$E-1) { $cont *= E_edge($degs->[$e], @$W) / $degs->[$e]; }
            for my $i (0..$n0-1) { my @dF; for my $e (0..$E-1) { push @dF, $degs->[$e] if $pairs[$e] % $n0 == $i; }
              $cont *= V_vertex(0, $gs->[$i], $W, @dF); last if $cont == 0; }
            next if $cont == 0;
            for my $j (0..$ninf-1) { my @dF; for my $e (0..$E-1) { push @dF, $degs->[$e] if int($pairs[$e] / $n0) == $j; }
              $cont *= V_vertex(1, $gs->[$n0 + $j], $W, @dF); last if $cont == 0; }
            # labelled weight: 1/(n0! ninf! E!) times number of degree orderings: we enumerate ordered compositions,
            # i.e. labelled edges with degrees, so weight 1/(n0! ninf! E!)
            $total += $cont / (fact($n0) * fact($ninf) * fact($E));
          }
        }
      } }
    }
  }
  return $total;
}
if (!caller) {
  my ($g, $d, @w) = @ARGV; my @W = map { Q($_) } @w;
  print "N_{$g,$d}(w = @w) = ", N_gd($g, $d, \@W), "\n";
}
1;
