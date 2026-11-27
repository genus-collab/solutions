# zeros_scan.pl h r Dmax : for each component k, list lattice points R (parts >= 2, sorted) with 2h-2 <= D(R) <= Dmax where Xi^(k)_h(R) = 0
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Comps;
my ($h, $r, $Dmax) = @ARGV; my $C = Comps::load("$FindBin::Bin/engine_h${h}_r${r}.txt");
my @pts; my $rec; $rec = sub { my ($pre, $min, $D, $k) = @_; if ($k == 0) { push @pts, [@$pre] if $D >= 2*$h - 2; return }
  for (my $j = $min; ($j - 1) + $D <= $Dmax; $j++) { $rec->([@$pre, $j], $j, $D + $j - 1, $k - 1); } }; $rec->([], 2, 0, $r);
for my $k (0..$h) { my @z; for my $p (@pts) { push @z, "{" . join(',', @$p) . "}" if Comps::comp_at($C, $h, $r, $k, $p) == 0; }
  print "genus $h r=$r k=$k: zeros with 2h-2 <= D <= $Dmax: @z\n"; }
# full Xi (all components zero)
