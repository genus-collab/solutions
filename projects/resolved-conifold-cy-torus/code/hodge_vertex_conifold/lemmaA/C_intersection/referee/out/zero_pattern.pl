use strict; use warnings; use lib '..'; use FindBin; use lib "$FindBin::Bin/.."; use rpoly; use Math::BigRat;
for my $hr ([2,1],[2,2],[3,1],[3,2],[4,1],[4,2]) { my ($h,$r) = @$hr; open(my $fh, '<', "$FindBin::Bin/ucomp_g${h}_r${r}.txt") or die; my $X0; while (<$fh>) { s/\r//; chomp; my @f = split /\s*\|\s*/; $X0 = parse_rp($f[3]) if $f[2] == 0; } close $fh;
  my %byD; if ($r == 1) { for my $j (1 .. 3*$h + 2) { my $D = $j - 1; push @{$byD{$D}}, [ "($j)", rp_eval($X0, $j) ]; } }
  else { for my $a (1 .. 3*$h+2) { for my $b ($a .. 3*$h+2) { my $D = $a + $b - 2; next if $D > 3*$h; push @{$byD{$D}}, [ "($a,$b)", rp_eval($X0, $a, $b) ]; } } }
  my @s; for my $D (sort { $a <=> $b } keys %byD) { my $nz = grep { !$_->[1]->is_zero } @{$byD{$D}}; push @s, "D=$D:" . ($nz ? "$nz/" . scalar(@{$byD{$D}}) . " nonzero" : "all 0"); }
  print "h=$h r=$r (parts >= 1): ", join('  ', @s), "\n"; }
