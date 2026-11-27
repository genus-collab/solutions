# rec2_test.pl hmax : Route B's (REC) for two leaves R = (j, y), tested as an identity of RATIONAL functions at random
# rational points (j, y) (not only lattice points), using the Hurwitz components (k = 0) of the engine tables:
#   (2h+1)(j+y-2h) Y_h(j,y) = Sum_{k=2}^{h+1} (2/(2k)!) (j+y-2(h+1-k))_(2k) Y_{h+1-k}(j,y)
#        + Sum_{h1+k1+h2+k2 = h+2} (2/(2k1)!)(2/(2k2)!) (j-2h1+1)_(2k1) (y-2h2+1)_(2k2) X_{h1}(j) X_{h2}(y),
# Y_g = Xi^(0)_g(j,y), X_g = Xi^(0)_g(j); genus 0: Y_0 = 1/(j+y), X_0(j) = 1/(j(j+1)).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my $hmax = shift // 4;
my $E = "$FindBin::Bin/nr_test";
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
sub ff { my ($x, $n) = @_; my $p = Q(1); $p *= ($x - $_) for 0..$n-1; return $p }
# r = 1 tables (k = 0)
my %X; for my $g (1..$hmax) { open my $fh, '<', "$FindBin::Bin/refdata/ref_engine_h${g}_r1.txt" or die; while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; next unless $k == 0; my @c;
    for my $t (split /\s*;\s*/, $rest) { my ($a, $v) = split /:/, $t; $c[$a] = Q($v); } $X{$g} = \@c; } close $fh; }
my %Y; for my $g (1..$hmax) { open my $fh, '<', "$E/engine_h${g}_r2.txt" or die "no r2 table $g"; while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; next unless $k == 0; my %c;
    for my $t (split /\s*;\s*/, $rest) { my ($m, $v) = split /:/, $t; $c{$m} = Q($v); } $Y{$g} = \%c; } close $fh; }
sub Xv { my ($g, $j) = @_; return Q(1) / ($j * ($j + 1)) if $g == 0; my $s = Q(0); my $c = $X{$g}; for my $a (0..$#$c) { next unless defined $c->[$a]; $s += $c->[$a] * $j ** $a; } return $s; }
sub Yv { my ($g, $j, $y) = @_; return Q(1) / ($j + $y) if $g == 0; my $s = Q(0); for my $m (keys %{$Y{$g}}) { my ($a, $b) = split /,/, $m; $s += $Y{$g}{$m} * $j ** $a * $y ** $b; } return $s; }
my @pts = ([Q('7/3'), Q('5/2')], [Q('-3/7'), Q('11/5')], [Q('13/4'), Q('-2/9')], [Q('9/2'), Q('7/2')], [Q('1/3'), Q('17/6')], [Q(3), Q(4)], [Q(5), Q(2)]);
for my $h (1..$hmax) { my $bad = 0;
  for my $pt (@pts) { my ($j, $y) = @$pt;
    my $lhs = (2*$h + 1) * ($j + $y - 2*$h) * Yv($h, $j, $y);
    my $rhs = Q(0);
    for my $k (2..$h+1) { my $hp = $h + 1 - $k; $rhs += Q(2) / fact(2*$k) * ff($j + $y - 2*$hp, 2*$k) * Yv($hp, $j, $y); }
    for my $h1 (0..$h+1) { for my $k1 (1..$h+1) { for my $h2 (0..$h+1) { my $k2 = $h + 2 - $h1 - $k1 - $h2; next if $k2 < 1;
          $rhs += Q(2) / fact(2*$k1) * Q(2) / fact(2*$k2) * ff($j - 2*$h1 + 1, 2*$k1) * ff($y - 2*$h2 + 1, 2*$k2) * Xv($h1, $j) * Xv($h2, $y); } } }
    $bad++ if $lhs != $rhs; }
  printf "genus %d: REC (two leaves) at %d rational points: %d failures\n", $h, scalar(@pts), $bad; }
