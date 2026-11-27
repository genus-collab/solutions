# t_tally.pl : collect the (A), (C), (M), (H) lines of all referee lattice outputs (lat*.txt, latm_*.txt, latf_*.txt) into one table
# per genus: Lemma-A points / violations, control points / non-zero, Hurwitz layers D = 2h-2, 2h-1, 2h, Conjecture-M lattice points.
use strict; use warnings; use FindBin;
my %T; my %src;
opendir(my $dh, $FindBin::Bin) or die; my @files = map { "$FindBin::Bin/$_" } sort grep { /^lat.*.txt$/ } readdir($dh); closedir $dh;
for my $f (@files) { open my $fh, '<', $f or next; while (<$fh>) {
    if (/^\(A\) genus (\d+) r = (\d+): (\d+) Lemma-A points .*?, (\d+) violations/) { $T{$1}{$2}{A} = [$3, $4]; }
    elsif (/^\(C\) genus (\d+) r = (\d+): (\d+) control points .*?, (\d+) non-zero; non-zero components: (.*)$/) { my ($h, $r, $n, $nz, $l) = ($1, $2, $3, $4, $5); my @c = split ' ', $l; my $all = grep { /:all$/ } @c; $T{$h}{$r}{C} = [$n, $nz, $all]; }
    elsif (/^\(M\) genus (\d+) r = (\d+): (\d+) lattice points .*?: (\d+) violations/) { $T{$1}{$2}{M} = [$3, $4]; }
    elsif (/^\(H\) genus (\d+) r = (\d+), D = (\d+): Hurwitz component zero at (\d+) points, non-zero at (\d+)/) { $T{$1}{$2}{"H$3"} = [$4, $5] if $4 + $5 > 0; } } close $fh; }
for my $h (sort { $a <=> $b } keys %T) { my ($nA, $bA, $nC, $nzC, $allC, $nM, $bM) = (0) x 7; my %H; my @rs = sort { $a <=> $b } keys %{$T{$h}};
  for my $r (@rs) { my $t = $T{$h}{$r}; if ($t->{A}) { $nA += $t->{A}[0]; $bA += $t->{A}[1]; } if ($t->{C}) { $nC += $t->{C}[0]; $nzC += $t->{C}[1]; $allC += $t->{C}[2]; }
    if ($t->{M}) { $nM += $t->{M}[0]; $bM += $t->{M}[1]; } for my $D (2*$h-2 .. 2*$h) { if ($t->{"H$D"}) { $H{$D}[0] += $t->{"H$D"}[0]; $H{$D}[1] += $t->{"H$D"}[1]; } } }
  printf "genus %d (r = %s): Lemma-A points %d, violations %d | controls %d, non-zero %d, all components non-zero %d | Conj-M lattice %d, violations %d | Hurwitz zero/non-zero: D=%d %d/%d, D=%d %d/%d, D=%d %d/%d\n",
    $h, join(',', @rs), $nA, $bA, $nC, $nzC, $allC, $nM, $bM, 2*$h-2, @{$H{2*$h-2} // [0,0]}, 2*$h-1, @{$H{2*$h-1} // [0,0]}, 2*$h, @{$H{2*$h} // [0,0]}; }
