package RefWK;
# Referee's independent implementation of the Witten-Kontsevich numbers <tau_{e_1} ... tau_{e_n}>_g (exact, Math::BigRat).
# Written from scratch (does not share code with ../WK.pm).  Algorithm:
#   string  <tau_0 tau_X>_g = Sum_j <tau_{X - e_j}>_g,   dilaton <tau_1 tau_X>_g = (2g-2+|X|) <tau_X>_g,
#   otherwise the DVV (Virasoro) recursion applied to the SMALLEST exponent k (>= 2) (the author applies it to the largest):
#   (2k+1)!! <tau_k tau_X>_g = Sum_{j in X} (2k+2e_j-1)!!/(2e_j-1)!! <tau_{k+e_j-1} tau_{X-j}>_g
#        + 1/2 Sum_{a+b=k-2} (2a+1)!! (2b+1)!! [ <tau_a tau_b tau_X>_{g-1} + Sum_{g1+g2=g, I+J=X} <tau_a tau_I>_{g1} <tau_b tau_J>_{g2} ].
# Base: <tau_0^3>_0 = 1, <tau_1>_1 = 1/24.  Unstable brackets are 0.
use strict; use warnings; use Math::BigRat; use Math::BigInt;
our %M;
my $Z = Math::BigRat->new(0);
my %DF;
sub df { my $m = shift; return Math::BigRat->new(1) if $m <= 1; return $DF{$m} //= do { my $r = Math::BigInt->new(1); for (my $k = $m; $k > 1; $k -= 2) { $r->bmul($k) } Math::BigRat->new($r) }; }
sub binom { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; my $r = 1; for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i } return $r; }
# sub-multisets of a sorted list: [ [I], [J], multiplicity ]
sub splits { my @x = @_; my (@v, @c); for my $y (@x) { if (@v && $v[-1] == $y) { $c[-1]++ } else { push @v, $y; push @c, 1 } }
  my @out = ([[], [], 1]);
  for my $i (0..$#v) { my @new; for my $o (@out) { for my $s (0..$c[$i]) {
        push @new, [[@{$o->[0]}, ($v[$i]) x $s], [@{$o->[1]}, ($v[$i]) x ($c[$i] - $s)], $o->[2] * binom($c[$i], $s)]; } } @out = @new; }
  return @out; }
sub wk {
  my ($g, @e) = @_;
  return $Z if $g < 0;
  my $n = scalar @e; return $Z if 2*$g - 2 + $n <= 0;
  my $s = 0; for (@e) { return $Z if $_ < 0; $s += $_ } return $Z if $s != 3*$g - 3 + $n;
  @e = sort { $a <=> $b } @e;
  my $key = "$g:" . join(',', @e);
  return $M{$key} if exists $M{$key};
  my $res;
  if ($g == 0 && $n == 3) { $res = Math::BigRat->new(1) }
  elsif ($g == 1 && $n == 1) { $res = Math::BigRat->new('1/24') }
  elsif ($e[0] == 0) { my @x = @e[1..$#e]; $res = Math::BigRat->new(0); my %seen;
    for my $j (0..$#x) { next if $x[$j] == 0; my @y = @x; $y[$j]--; $res += wk($g, @y); } }
  elsif ($e[0] == 1) { my @x = @e[1..$#e]; $res = wk($g, @x) * (2*$g - 2 + scalar(@x)); }
  else { my $k = $e[0]; my @x = @e[1..$#e]; $res = Math::BigRat->new(0);
    for my $j (0..$#x) { my @y = @x; my $ej = $y[$j]; $y[$j] = $k + $ej - 1;
      $res += df(2*$k + 2*$ej - 1) / df(2*$ej - 1) * wk($g, @y); }
    my $half = Math::BigRat->new('1/2');
    for my $a (0..$k-2) { my $b = $k - 2 - $a; my $cf = df(2*$a + 1) * df(2*$b + 1); my $t = wk($g - 1, $a, $b, @x);
      for my $sp (splits(@x)) { my ($I, $J, $mult) = @$sp;
        for my $g1 (0..$g) { my $v1 = wk($g1, $a, @$I); next if $v1 == 0; my $v2 = wk($g - $g1, $b, @$J); next if $v2 == 0; $t += $v1 * $v2 * $mult; } }
      $res += $half * $cf * $t; }
    $res = $res / df(2*$k + 1); }
  $M{$key} = $res; return $res; }
1;
