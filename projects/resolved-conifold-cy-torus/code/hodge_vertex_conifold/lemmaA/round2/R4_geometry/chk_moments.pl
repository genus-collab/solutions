# chk_moments.pl h... : reproduce the component tables refdata/ref_engine_h<h>_r1.txt from the moments (Pair.pm).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Pair;
sub Q { Math::BigRat->new(shift) }
for my $h (@ARGV) {
  my @comp = Pair::omega_components($h);
  my @T; open my $fh, '<', "$FindBin::Bin/refdata/ref_engine_h${h}_r1.txt" or die; while (my $l = <$fh>) { chomp $l; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $l, 4; $rest //= ''; my @c;
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($a, $v) = split /:/, $t; $c[$a] = Q($v); } $T[$k] = \@c; } close $fh;
  my $bad = 0;
  for my $k (0..$h) { my $p = Pair::jpoly($h, $comp[$k]); my $n = $#$p > $#{$T[$k]} ? $#$p : $#{$T[$k]};
    for my $a (0..$n) { my $x = $p->[$a] // Q(0); my $y = $T[$k][$a] // Q(0); $bad++ if $x != $y; } }
  print "genus $h: components from moments vs table: ", ($bad ? "$bad MISMATCHES" : "all equal"), "\n"; }
