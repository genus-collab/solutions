# hurwitz_extra.pl h : the extra layer of zeros of the Hurwitz component (k = 0, the single class Lambda^vee(1), i.e. c(T(-log))):
#   Xi^(0)_h(R) = 0 for every multiset R of integers >= 2 with D(R) = 2h-1 ?   (r = 1: the extra zero j = 2h)
# also reports the other components there (expected non-zero) and the Hurwitz component at D = 2h-2 and 2h (expected non-zero).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use Time::HiRes qw(time);
$| = 1;
my $h = shift; my $t0 = time; Hodge::setup($h);
my @pts; my $rec; $rec = sub { my ($pre, $min, $D) = @_; push @pts, [@$pre] if @$pre && $D >= 2*$h - 2;
  for (my $j = $min; ($j - 1) + $D <= 2*$h; $j++) { $rec->([@$pre, $j], $j, $D + $j - 1); } };
$rec->([], 2, 0);
my %cnt;
for my $p (sort { scalar(@$a) <=> scalar(@$b) || join(',', @$a) cmp join(',', @$b) } @pts) { my $D = 0; $D += $_ - 1 for @$p;
  my $v = HodgePt::WP([], $p); my $z0 = (($v->[0] // 0) == 0) ? 1 : 0; my @nzk = grep { ($v->[$_] // 0) != 0 } 0..$h;
  $cnt{"$D|$z0"}++;
  printf "genus %d R = {%s} D = %d: Hurwitz component %s; non-zero components {%s}\n", $h, join(',', @$p), $D, ($z0 ? "ZERO" : "non-zero"), join(' ', @nzk); }
for my $D (2*$h-2 .. 2*$h) { printf "SUMMARY genus %d, D = %d: %d points with Xi^(0) = 0, %d with Xi^(0) != 0\n", $h, $D, $cnt{"$D|1"} // 0, $cnt{"$D|0"} // 0; }
printf "time %.1f s\n", time - $t0;
