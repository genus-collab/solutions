# loc_P1.pl -- GP localization of Int_{[Mbar_g(P^1,d)]^vir} 1 (= 0 when vdim = 2g-2+2d > 0) with weights (1, -1),
# and of the tau -> oo limit of the X_3 graph contributions (report.md, Proposition D3.3):
#   Cont_{X3}(Gamma) = (-1)^{g-1+d} tau^{2(g-1+d)} Cont_{P1}(Gamma) (1 + O(1/tau)).
# Usage: perl loc_P1.pl g d
use strict; use warnings; use File::Basename qw(dirname); { my $dd = dirname(__FILE__); require "$dd/loc_conifold.pl"; }
no warnings 'redefine';
our $MODE = "P1";
my $w = Q(1);
sub E_P1 { my $d = shift; return Q($d)**(2*$d) / (Q((-1)**$d) * fact($d)**2 * $w**(2*$d)); }
my %VP;
sub V_P1 { my ($p, $g, @dF) = @_; my $t = $p == 0 ? $w : -$w; my $val = scalar @dF;
  my $key = join("|", $p, $g, sort { $a <=> $b } @dF); return $VP{$key} if exists $VP{$key};
  my $res;
  if ($g == 0 && $val == 1) { $res = $t / $dF[0]; }
  elsif ($g == 0 && $val == 2) { $res = $t / ($t / $dF[0] + $t / $dF[1]); }
  else { my $cls = omega_class($g, $t); my $N = 3*$g - 3 + $val;
    my @ser; for my $dd (@dF) { my @s; my $c = Q($dd) / $t; for my $a (0..$N) { push @s, $c * (Q($dd) / $t)**$a; } push @ser, \@s; }
    $res = $t**($val - 1) * $t**$g * vertex($g, $cls, @ser); }
  $VP{$key} = $res; return $res; }
sub N_P1 { my ($g, $d) = @_; my $total = Q(0);
  for my $E (1..$d) { for my $degs (compositions($d, $E)) { for my $n0 (1..$E) { for my $ninf (1..$E) {
    next if $n0 + $ninf > $E + 1; my $h1 = $E - $n0 - $ninf + 1; next if $h1 < 0 || $h1 > $g; my $np = $n0 * $ninf;
    for my $code (0..$np**$E - 1) { my @pairs; my $c = $code; for my $e (0..$E-1) { push @pairs, $c % $np; $c = int($c / $np); }
      my @v0 = (0) x $n0; my @vI = (0) x $ninf; my @par = (0..$n0+$ninf-1);
      my $find = sub { my $x = shift; while ($par[$x] != $x) { $x = $par[$x]; } return $x; };
      for my $e (0..$E-1) { my $i = $pairs[$e] % $n0; my $j = int($pairs[$e] / $n0); $v0[$i]++; $vI[$j]++; my ($a, $b) = ($find->($i), $find->($n0 + $j)); $par[$a] = $b if $a != $b; }
      next if grep { $_ == 0 } (@v0, @vI); my $root = $find->(0); next if grep { $find->($_) != $root } 0..$n0+$ninf-1;
      for my $gs (genus_assign($g - $h1, $n0 + $ninf)) { my $cont = Q(1);
        for my $e (0..$E-1) { $cont *= E_P1($degs->[$e]) / $degs->[$e]; }
        for my $i (0..$n0-1) { my @dF = map { $degs->[$_] } grep { $pairs[$_] % $n0 == $i } 0..$E-1; $cont *= V_P1(0, $gs->[$i], @dF); }
        for my $j (0..$ninf-1) { my @dF = map { $degs->[$_] } grep { int($pairs[$_] / $n0) == $j } 0..$E-1; $cont *= V_P1(1, $gs->[$n0 + $j], @dF); }
        $total += $cont / (fact($n0) * fact($ninf) * fact($E)); } } } } } }
  return $total; }
my ($g, $d) = @ARGV;
print "Int_[Mbar_$g(P^1,$d)]^vir 1 by localization = ", N_P1($g, $d), "\n";
# large-tau scaling of the X_3 sum: N_{g,d}(tau) is constant, but check graph-by-graph scaling via total of |top coefficient|
