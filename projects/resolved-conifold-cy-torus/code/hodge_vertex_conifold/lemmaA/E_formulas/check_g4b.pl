# check_g4b.pl : engine (HodgePt) vs the Mathematica values of hvc_multiprobe_g4b_out.txt (genus 4, r = 3..5), at x = 2 and x = 1/3
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use MParse; use Poly;
sub R { Math::BigRat->new(shift) }
Hodge::setup(4);
open my $fh, '<', "$FindBin::Bin/../../hvc_multiprobe_g4b_out.txt" or die; my ($n, $ok) = (0, 0);
while (<$fh>) { next unless /^genus 4, R = \{([\d, ]+)\}, D = \d+ \([^)]*\): (.*?)\s+\[/; my @R = split /,\s*/, $1; my $val = $2; my $r = scalar @R;
  my $v = HodgePt::WP([], \@R);
  for my $x (R(2), R("1/3")) { my $tau = -$x / (1 + $x); my $sig = $tau * (1 + $tau); my $vv = R(-1)/$sig; my $s = R(0); for my $k (0..$#$v) { $s += ($v->[$k] // 0) * $vv ** $k; }
    my $mine = $s / (1 + $x) ** (9 + 2*$r);
    my $theirs = $val eq '0' ? R(0) : do { (my $e = $val) =~ s/^non-zero, //; Poly::const_val(MParse::parse($e, { x => $x })) };
    $n++; $ok++ if $mine == $theirs; print "R = {@R}, x = $x: ", ($mine == $theirs ? "agree" : "DISAGREE"), "\n"; } }
print "hvc_multiprobe_g4b: $ok / $n comparisons agree\n";
