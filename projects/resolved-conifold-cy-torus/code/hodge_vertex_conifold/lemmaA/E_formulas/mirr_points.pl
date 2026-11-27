# mirr_points.pl h : exact check of M_irr(h) (Lemma A at genus h) at every irreducible point, by point evaluation (HodgePt.pm):
# every multiset R of integers >= 2 with D(R) <= 2h-3 (all r <= 2h-3), all components Xi^(k), k = 0..h;
# plus the control points D(R) = 2h-2 (expected non-zero).  Prints one line per point.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use Time::HiRes qw(time);
$| = 1; my $h = shift; my $t0 = time; Hodge::setup($h); printf "T-table setup done %.1f s
", time - $t0;
my @pts; my $rec; $rec = sub { my ($pre, $min, $D) = @_; push @pts, [@$pre] if @$pre;
  for (my $j = $min; ($j - 1) + $D <= 2*$h - 2; $j++) { $rec->([@$pre, $j], $j, $D + $j - 1); } };
$rec->([], 2, 0);
my ($nA, $bad, $nC, $nz) = (0, 0, 0, 0);
for my $p (sort { scalar(@$a) <=> scalar(@$b) || join(',', @$a) cmp join(',', @$b) } @pts) { my $D = 0; $D += $_ - 1 for @$p;
  my $v = HodgePt::WP([], $p); my @nzk = grep { ($v->[$_] // 0) != 0 } 0..$h;
  if ($D <= 2*$h - 3) { $nA++; $bad++ if @nzk; printf "genus %d R = {%s} D = %d: %s\n", $h, join(',', @$p), $D, (@nzk ? "NON-ZERO components @nzk" : "all components 0"); }
  else { $nC++; $nz++ if @nzk; printf "genus %d R = {%s} D = %d (control): non-zero components {%s}\n", $h, join(',', @$p), $D, join(' ', @nzk); } }
printf "SUMMARY genus %d: %d irreducible points, %d with a non-zero component; %d control points (D = 2h-2), %d non-zero; %.1f s\n", $h, $nA, $bad, $nC, $nz, time - $t0;
