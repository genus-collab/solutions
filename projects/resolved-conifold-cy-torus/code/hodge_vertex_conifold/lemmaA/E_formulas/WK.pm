package WK;
# Witten-Kontsevich intersection numbers <tau_{a_1} ... tau_{a_n}>_g on Mbar_{g,n}, exact (Math::BigRat).
# String and dilaton equations first, then the DVV (Virasoro) recursion
#   <tau_{k+1} prod_X tau_d>_g = 1/(2k+3)!! [ Sum_j (2k+2d_j+1)!!/(2d_j-1)!! <tau_{d_j+k} prod_{X-j}>_g
#        + 1/2 Sum_{r+s=k-1} (2r+1)!!(2s+1)!! ( <tau_r tau_s prod_X>_{g-1} + Sum_{g1+g2=g, I+J=X} <tau_r prod_I>_{g1} <tau_s prod_J>_{g2} ) ].
# Base values <tau_0^3>_0 = 1, <tau_1>_1 = 1/24.
use strict; use warnings;
use Math::BigRat;
our %memo;
sub R { Math::BigRat->new(shift) }
my $ZERO = R(0);
sub dfact { # (2m-1)!! for m >= 0, as BigRat ; dfact(m) = (2m-1)!!
  my $m = shift; my $r = R(1); for (my $k = 2*$m-1; $k > 1; $k -= 2) { $r *= $k } return $r; }
my %dfc;
sub df { my $m = shift; $dfc{$m} //= dfact($m); return $dfc{$m}; }

sub wk {
  my ($g, @a) = @_;
  return $ZERO if $g < 0;
  for (@a) { return $ZERO if $_ < 0 }
  my $n = scalar @a;
  return $ZERO if 2*$g - 2 + $n <= 0;
  my $s = 0; $s += $_ for @a;
  return $ZERO if $s != 3*$g - 3 + $n;
  @a = sort { $a <=> $b } @a;
  my $key = "$g|" . join(',', @a);
  return $memo{$key} if exists $memo{$key};
  my $res;
  if ($g == 0 && $n == 3) { $res = R(1); }
  elsif ($g == 1 && $n == 1) { $res = R("1/24"); }
  elsif ($a[0] == 0) {        # string
    my @b = @a[1..$#a]; $res = R(0);
    for my $j (0..$#b) { next if $b[$j] == 0; my @c = @b; $c[$j]--; $res += wk($g, @c); }
  }
  elsif ($a[0] == 1) {        # dilaton
    my @b = @a[1..$#a];
    $res = wk($g, @b) * (2*$g - 2 + scalar(@b));
  }
  else {                      # DVV on the largest entry
    my $kp1 = pop @a; my $k = $kp1 - 1; my @X = @a;
    $res = R(0);
    for my $j (0..$#X) {
      my @c = @X; my $d = $c[$j]; $c[$j] = $d + $k;
      my $coef = df($k + $d + 1) / df($d);     # (2k+2d+1)!! / (2d-1)!!
      $res += $coef * wk($g, @c);
    }
    my $half = R("1/2");
    for my $r (0..$k-1) { my $sx = $k - 1 - $r;
      my $cf = df($r+1) * df($sx+1);           # (2r+1)!! (2s+1)!!
      my $t = wk($g-1, $r, $sx, @X);
      # splittings over sub-multisets of X
      my @subs = submultisets(@X);
      for my $sb (@subs) {
        my ($I, $J, $mult) = @$sb;
        for my $g1 (0..$g) { my $g2 = $g - $g1;
          my $v1 = wk($g1, $r, @$I); next if $v1 == 0;
          my $v2 = wk($g2, $sx, @$J); next if $v2 == 0;
          $t += $mult * $v1 * $v2;
        }
      }
      $res += $half * $cf * $t;
    }
    $res = $res / df($k + 2);                  # (2k+3)!!
  }
  $memo{$key} = $res;
  return $res;
}

# all splittings of the multiset X into (I, J) with multiplicity prod binom(m_v, i_v)
sub submultisets {
  my @X = @_;
  my %t; $t{$_}++ for @X;
  my @vals = sort { $a <=> $b } keys %t;
  my @out = ([[], [], 1]);
  for my $v (@vals) {
    my $m = $t{$v}; my @new;
    for my $o (@out) {
      for my $i (0..$m) {
        my @I = (@{$o->[0]}, ($v) x $i); my @J = (@{$o->[1]}, ($v) x ($m - $i));
        push @new, [\@I, \@J, $o->[2] * binom($m, $i)];
      }
    }
    @out = @new;
  }
  return @out;
}
sub binom { my ($n, $k) = @_; return 0 if $k < 0 || $k > $n; my $r = 1; for my $i (1..$k) { $r = $r * ($n - $k + $i) / $i } return $r; }

1;
