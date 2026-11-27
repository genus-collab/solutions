# mpoly.pl : collect the slice series C^{(m)}_s (slice_m<m>.txt, m = 1..11; m = 0: C_0 = S1S2) and test whether each coefficient
# [lambda^{2n} sigma^p] C^{(m)}_s is a polynomial in m (exact Newton interpolation; the degree must be < number of points - 1 to be a test).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use SerS;
sub R { Math::BigRat->new(shift) }
my %C;   # C{m}{s}{n}{p}
for my $m (1..16) { open my $fh, '<', "$FindBin::Bin/slice_m$m.txt" or next;
  while (<$fh>) { next unless /C_(\d+) \[lambda\^(\d+)\] = (.*)$/; my ($s, $n, $rest) = ($1, $2/2, $3); next if $rest =~ /^0\s*$/;
    while ($rest =~ /\(([-\d\/]+)\)s\^(\d+)/g) { $C{$m}{$s}{$n}{$2} = R($1); } } }
# m = 0: S1S2
my $S12 = SerS::S1S2(); for my $n (0..5) { my $c = SerS::coeff($S12, $n); for my $p (0..$#$c) { $C{0}{0}{$n}{$p} = $c->[$p] if $c->[$p] != 0; } }
my @ms = sort { $a <=> $b } keys %C;
print "slices available: m = @ms\n";
for my $s (0..10) { for my $n (0..5) { for my $p (0..$n) {
  my @x; my @y; for my $m (@ms) { push @x, R($m); push @y, ($C{$m}{$s}{$n}{$p} // R(0)); }
  my $allzero = 1; $allzero = 0 if grep { $_ != 0 } @y; next if $allzero;
  # Newton divided differences
  my @dd = @y; my $N = scalar @x; for my $k (1..$N-1) { for (my $i = $N-1; $i >= $k; $i--) { $dd[$i] = ($dd[$i] - $dd[$i-1]) / ($x[$i] - $x[$i-$k]); } }
  my $deg = 0; for my $k (0..$N-1) { $deg = $k if $dd[$k] != 0; }
  my $test = $deg < $N - 1 ? "polynomial of degree $deg (tested with $N points)" : "degree >= " . ($N-1) . " (NO TEST)";
  printf "  [lambda^%d sigma^%d] C_%d(m): %s\n", 2*$n, $p, $s, $test; } } }
