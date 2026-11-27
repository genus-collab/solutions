use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Comps; use Math::BigRat;
my ($h, $r, @pts) = @ARGV; my $C = Comps::load("$FindBin::Bin/engine_h${h}_r${r}.txt");
for my $p (@pts) { my @j = split /,/, $p; print "Xi^(k)_$h(@j): ", join("  ", map { "k=$_: " . Comps::comp_at($C, $h, $r, $_, \@j) } 0..$h), "\n"; }
