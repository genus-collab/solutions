# cyH_test.pl gmax : does the boundary-sum identity H(g) of the lambda_1 proof hold for the FULL CY class, component by component?
#   B^(k)_g(j) := [y^1] Xi^(k)_g(j,y)          (two-leaf tables nr_test/engine_h<g>_r2.txt, component k)
#   BS^(k)_g(j) := Sum_{g1=0}^{g-1} Sum_{k1+k2=k} Xi^(k2)_{g-g1}(0) Y^(k1)_{g1}(j),  Y^(k)_{g1}(j) = (j-2g1+1) Xi^(k)_{g1}(j)  (g1 >= 1),
#                  Y^(0)_0 = 1/j, Y^(k>0)_0 = 0  (the separating boundary sum; Omega is multiplicative on separating divisors).
# Prints the integers j in 1..3g where B^(k)_g(j) = BS^(k)_g(j).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my $gmax = shift // 5;
sub load1 { my $g = shift; my @T; open my $fh, '<', "$FindBin::Bin/refdata/ref_engine_h${g}_r1.txt" or die; while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; $rest //= ''; my @c;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($a, $v) = split /:/, $t; $c[$a] = Q($v); } $T[$k] = \@c; } close $fh; return \@T; }
sub load2y1 { my $g = shift; my @T; open my $fh, '<', "$FindBin::Bin/nr_test/engine_h${g}_r2.txt" or die; while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; $rest //= ''; my @c;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $v) = split /:/, $t; my ($a, $b) = split /,/, $m; $c[$a] = Q($v) if $b == 1; } $T[$k] = \@c; } close $fh; return \@T; }
sub ev { my ($p, $x) = @_; my $s = Q(0); return $s unless $p; for my $a (0..$#$p) { next unless defined $p->[$a]; $s += $p->[$a] * Q($x) ** $a; } return $s; }
my %X; $X{$_} = load1($_) for 1..$gmax;
for my $g (2..$gmax) { my $B = load2y1($g);
  for my $k (0..$g) { my @eq; my @ne;
    for my $j (1..3*$g) { my $bs = Q(0);
      for my $g1 (0..$g-1) { my $g2 = $g - $g1;
        for my $k1 (0..$k) { my $k2 = $k - $k1; next if $k2 > $g2; my $c2 = ev($X{$g2}[$k2], 0); next if $c2 == 0;
          my $y; if ($g1 == 0) { $y = $k1 == 0 ? Q(1)/$j : Q(0); } else { next if $k1 > $g1; $y = ($j - 2*$g1 + 1) * ev($X{$g1}[$k1], $j); }
          $bs += $c2 * $y; } }
      if (ev($B->[$k], $j) == $bs) { push @eq, $j } else { push @ne, $j } }
    printf "g=%d k=%d: B = BS at j = %s\n", $g, $k, (@eq ? join(',', @eq) : 'none'); } }
