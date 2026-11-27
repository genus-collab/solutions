# t_g4b.pl : the referee engine against the Mathematica genus-4 multi-leaf values of ../../../hvc_multiprobe_g4b_out.txt (r = 3..5),
# symbolic in x there; compared here at x = 2, 1/3, -3/5, 7 (c-normalisation: Xi^c = c^{-(3h-3+2r)} Xi^{c=1}, tau = -x/(1+x)).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge; use RefParse;
sub Q { Math::BigRat->new(shift) }
my $f = "$FindBin::Bin/../../../hvc_multiprobe_g4b_out.txt"; open my $fh, '<', $f or die "no $f"; my @lines = grep { /^genus 4, R/ } <$fh>; close $fh;
my $h = 4; my ($n, $bad) = (0, 0);
for my $x (map { Q($_) } qw(2 1/3 -3/5 7)) { my $c = 1 + $x; my $tau = -$x / $c; my $om = RefHodge::omega($h, [Q(1), $tau, -1 - $tau], 3*$h - 3);
  for my $line (@lines) { my ($Rs, $val) = $line =~ /R = \{([\d, ]+)\}.*?\):\s+(.*?)\s+\[/ or next; my @R = split /,\s*/, $Rs; my $r = scalar @R;
    my $mine = RefHodge::xi_value(RefHodge::xi($h, $om, \@R)) / ($c ** (3*$h - 3 + 2*$r));
    my $theirs = $val eq '0' ? Q(0) : do { (my $e = $val) =~ s/^non-zero,\s*//; RefParse::evalx($e, { x => $x }) };
    $n++; $bad++ if $mine != $theirs; } }
printf "genus 4 multi-leaf Mathematica values (13 points x 4 values of x): %d comparisons, %d mismatches\n", $n, $bad;
