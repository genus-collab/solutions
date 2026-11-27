# ref_n3_k2.pl : N3 (CY analogue of M3) and the K2 reduction of the k = 2 component, re-checked on the MATHEMATICA component
# tables extracted by data_components.pl (genus 2, 3, and 4 when available); genus-1 one-leaf components and the genus h-1 three-leaf
# Hurwitz tables (needed only for the irreducible boundary of ch_3) are taken from Route E's engine tables.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my $E = "$FindBin::Bin/../../../E_formulas"; my $O = "$FindBin::Bin/out";
my @H = grep { -s "$O/data_h${_}_r2.txt" } (2, 3, 4);
sub load { my $f = shift; open my $fh, '<', $f or die "no $f"; my @C; while (my $l = <$fh>) { chomp $l; next unless $l =~ /\|/; my ($h, $r, $k, $rest) = split /\s*\|\s*/, $l, 4; my %c; $rest //= '';
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $v) = split /:/, $t; $m =~ s/\s//g; $c{$m} = Q($v); } $C[$k] = \%c; } close $fh; return \@C; }
my (%R1, %R2, %R3);
$R1{1} = load("$E/engine_h1_r1.txt"); $R3{1} = load("$E/engine_h1_r3.txt"); $R3{2} = load("$E/engine_h2_r3.txt"); $R3{3} = load("$E/engine_h3_r3.txt");
for my $h (2, 3, 4) { $R1{$h} = load("$O/data_h${h}_r1.txt") }
for my $h (@H) { $R2{$h} = load("$O/data_h${h}_r2.txt") }
$R2{1} = load("$E/engine_h1_r2.txt");
sub e1 { my ($g, $k, $j) = @_; my $c = $R1{$g}[$k] or return Q(0); my $s = Q(0); for my $m (keys %$c) { $s += $c->{$m} * Q($j) ** $m } return $s }
sub e1c { my ($g, $k, $b) = @_; my $c = $R1{$g}[$k] or return Q(0); return $c->{$b} // Q(0) }
sub e2y { my ($g, $k, $j, $cc) = @_; my $C = $R2{$g}[$k] or return Q(0); my $s = Q(0); for my $m (keys %$C) { my ($a, $b) = split /,/, $m; next unless $b == $cc; $s += $C->{$m} * Q($j) ** $a } return $s }
sub e3yz { my ($g, $j, $cy, $cz) = @_; my $C = $R3{$g}[0] or return Q(0); my $s = Q(0); for my $m (keys %$C) { my ($a, $b, $d) = split /,/, $m; next unless $b == $cy && $d == $cz; $s += $C->{$m} * Q($j) ** $a } return $s }
# ---- N3 ----
for my $g (@H) { for my $k (0..$g) { my (@eq, @ne);
    for my $j (1..3*$g) { my $bs = Q(0);
      for my $g1 (0..$g-1) { my $g2 = $g - $g1; for my $k1 (0..$k) { my $k2 = $k - $k1; next if $k2 > $g2; my $c2 = e1($g2, $k2, 0); next if $c2 == 0;
          my $y; if ($g1 == 0) { $y = $k1 == 0 ? Q(1)/$j : Q(0) } else { next if $k1 > $g1; $y = ($j - 2*$g1 + 1) * e1($g1, $k1, $j) } $bs += $c2 * $y } }
      if (e2y($g, $k, $j, 1) == $bs) { push @eq, $j } else { push @ne, $j } }
    printf "N3 (Mathematica data) g=%d k=%d: [y]Xi^(k)(j,y) = separating sum at j = %s\n", $g, $k, (@eq ? join(',', @eq) : 'none'); } }
# ---- K2 ----
for my $h (@H) { my ($n, $bad) = (0, 0);
  for my $j (1..3*$h) { my $jq = Q($j);
    my $A = (2*$h - 1) * e1($h, 1, $j) - e2y($h, 1, $j, 1) - (e1($h, 1, $j) - e1($h, 1, 0)) / $jq;
    $A += Q(1)/2 * ($jq - 2*$h + 2) * ($jq - 2*$h + 3) * e1($h - 1, 1, $j);
    for my $g1 (1..$h-1) { my $g2 = $h - $g1; $A += ($jq - 2*$g1 + 1) * (e1($g1, 1, $j) * e1($g2, 0, 0) + e1($g1, 0, $j) * e1($g2, 1, 0)) }
    my $k3 = (2*$h - 1) * e1($h, 0, $j) - e2y($h, 0, $j, 1) - e2y($h, 0, $j, 2) - e2y($h, 0, $j, 3);
    my $p3 = Q(0); for my $a (3..3*$h) { $p3 += e1c($h, 0, $a) * $jq ** ($a - 3) }
    my $irr = Q(1)/2 * (2 * e3yz($h - 1, $j, 2, 0) - e3yz($h - 1, $j, 1, 1));
    my $sep = Q(0); for my $g1 (1..$h-1) { my $g2 = $h - $g1; my $Y0 = ($jq - 2*$g1 + 1) * e1($g1, 0, $j);
      die "need two-leaf genus-$g1 data" unless exists $R2{$g1}; my ($Y1, $Y2) = (e2y($g1, 0, $j, 1), e2y($g1, 0, $j, 2));
      $sep += $Y2 * e1c($g2, 0, 0) - $Y1 * e1c($g2, 0, 1) + $Y0 * e1c($g2, 0, 2) }
    my $pred = $A / 24 + ($k3 - $p3 + $irr + $sep) / 120; $n++; $bad++ if $pred != e1($h, 2, $j) }
  printf "K2 reduction (Mathematica one/two-leaf data, Route E three-leaf genus %d): genus %d, %d values of j: %s\n", $h - 1, $h, $n, ($bad ? "$bad MISMATCHES" : "all equal"); }
