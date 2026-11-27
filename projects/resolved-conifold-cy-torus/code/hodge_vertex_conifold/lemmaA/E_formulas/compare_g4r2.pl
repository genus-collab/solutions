# compare_g4r2.pl : compare the Mathematica output of E_crosscheck_g4r2.wls (symbolic Xi[4,2], Xi[3,3], c-normalisation)
# with the Perl engine files engine_h4_r2.txt and engine_h3_r3.txt (components Xi^(k), c = 1).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Poly; use MParse; use Comp; use Comps;
open my $fh, '<', "$FindBin::Bin/E_crosscheck_g4r2_out.txt" or die "run E_crosscheck_g4r2.wls first";
while (<$fh>) { next unless /^genus (\d) r=(\d) symbolic: (.*)$/; my ($h, $r, $e) = ($1, $2, $3);
  my ($comp, $ok) = Comp::components($e, $h, $r);
  my $C = Comps::load("$FindBin::Bin/engine_h${h}_r${r}.txt"); my ($n, $bad) = (0, 0);
  for my $k (0..$h) { my $p = $comp->[$k]; my %seen;
    for my $m (keys %$p) { my $ex = Poly::mono_parse($m); my $key = join(',', map { $ex->{"jj$_"} // 0 } 1..$r); $seen{$key} = 1; $n++;
      $bad++ if ($C->{$h}{$r}{$k}{$key} // Math::BigRat->new(0)) != $p->{$m}; }
    for my $key (keys %{$C->{$h}{$r}{$k} // {}}) { next if $seen{$key}; $n++; $bad++ if $C->{$h}{$r}{$k}{$key} != 0; } }
  print "genus $h r=$r: 1/sigma-structure ", ($ok ? "ok" : "FAILED"), "; $n coefficients compared with the engine, $bad mismatches\n"; }
