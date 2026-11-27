# t_slicepoly.pl : check the author's claim "slice series C_s(m) are polynomials in m divisible by C(m, s+1) (17 slices, degrees <= 10)"
# directly from the stored slice files ../slice_m0..16.txt (m = 0 and 1: C_0 = S1 S2, C_s = 0 otherwise; read from the files where present).
# For every (s, lambda-order n, sigma-power p): the 17 values y(m), m = 0..16; degree of the interpolating polynomial (via finite differences)
# and vanishing at m = 0..s.
use strict; use warnings; use FindBin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my %C; my %seen;
for my $m (0..16) { my $f = "$FindBin::Bin/../slice_m$m.txt"; open my $fh, '<', $f or next;
  while (<$fh>) { next unless /C_(\d+) \[lambda\^(\d+)\] = (.*)$/; my ($s, $n, $rest) = ($1, $2/2, $3); $seen{"$s|$n"} = 1; while ($rest =~ /\(([-\d\/]+)\)s\^(\d+)/g) { $C{$m}{$s}{$n}{$2} = Q($1); } } }
print "slice m = 0 file content: ", (-e "$FindBin::Bin/../slice_m0.txt" ? `cat "$FindBin::Bin/../slice_m0.txt"` : "missing"), "\n";
# m = 0 (dilaton slice): C_0 = S1 S2 exactly (string/dilaton), C_s = 0 for s >= 1 (as in ../HFun.pm); computed here
{ use lib "$FindBin::Bin"; require RefHodge; my $K = 6; my @c = (Q(0)); for my $n (1..$K) { my $f = Q(1); $f *= $_ for 2..2*$n; push @c, -RefHodge::bern(2*$n) / (2*$n * $f); }
  my @pp = ([Q(2)], [Q(1)]); my $padd = sub { my ($a, $b) = @_; my $nn = $#$a > $#$b ? $#$a : $#$b; [ map { ($a->[$_] // Q(0)) + ($b->[$_] // Q(0)) } 0..$nn ] };
  my $pmul = sub { my ($a, $b) = @_; my @r = map { Q(0) } 0..($#$a + $#$b); for my $i (0..$#$a) { for my $j (0..$#$b) { $r[$i+$j] += $a->[$i] * $b->[$j] } } \@r };
  for my $n (2..2*$K) { push @pp, $padd->($pp[$n-1], $pmul->([Q(0), Q(1)], $pp[$n-2])); }
  my @lg = ([Q(0)], map { [ map { $_ * $c[$_[0] // 0] } () ] } ()); @lg = ([Q(0)]); for my $n (1..$K) { push @lg, [ map { $_ * $c[$n] } @{$pp[2*$n]} ]; }
  my @ex = ([Q(1)]); for my $n (1..$K) { my $s = [Q(0)]; for my $mm (1..$n) { my $t = $pmul->($lg[$mm], $ex[$n-$mm]); $s = $padd->($s, [ map { $_ * $mm } @$t ]); } $ex[$n] = [ map { $_ / $n } @$s ]; }
  for my $n (0..$K) { for my $p (0..$#{$ex[$n]}) { $C{0}{0}{$n}{$p} = $ex[$n][$p] if $ex[$n][$p] != 0; } } }
my ($nchk, $degmax, $divbad, $nonpoly) = (0, 0, 0, 0); my %maxdeg;
for my $key (sort keys %seen) { my ($s, $n) = split /\|/, $key; next if $n == 0;
  for my $p (0..$n+1) { my @y = map { $C{$_}{$s}{$n}{$p} // Q(0) } 0..16; next unless grep { $_ != 0 } @y;
    # finite differences
    my @d = @y; my $deg = -1; for my $k (0..16) { if (grep { $_ != 0 } @d) { $deg = $k } else { last } my @e; for my $i (0..$#d-1) { push @e, $d[$i+1] - $d[$i]; } @d = @e; last unless @d; }
    # deg = largest k with Delta^k y not identically zero on the available points (if deg <= 15 the data determine a polynomial of that degree)
    $nchk++; $degmax = $deg if $deg > $degmax; $maxdeg{$s} = $deg if !defined $maxdeg{$s} || $deg > $maxdeg{$s}; $nonpoly++ if $deg > 10;
    for my $m (0..$s) { if ($y[$m] != 0 && $m >= 2) { $divbad++; print "  C_$s [lambda^", 2*$n, "] sigma^$p non-zero at m = $m\n"; } } } }
printf "checked %d (s, order, sigma-power) coefficient sequences; maximal degree in m = %d (claimed <= 10); %d sequences of degree > 10; %d violations of vanishing at 2 <= m <= s\n", $nchk, $degmax, $nonpoly, $divbad;
print "max degree per s: ", join(' ', map { "s=$_:$maxdeg{$_}" } sort { $a <=> $b } keys %maxdeg), "\n";
