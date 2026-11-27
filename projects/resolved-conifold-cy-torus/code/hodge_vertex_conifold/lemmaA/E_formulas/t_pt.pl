use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt; use Comps;
for my $case ([3, [2,3]], [3, [0,5]], [4, [2,3,4]], [4, [1,2,5]], [4, [2,2,2,3,5]], [5, [3,7]], [2, [0,1,4,4]]) {
  my ($h, $J) = @$case; Hodge::setup($h); %HodgePt::memo = (); %Hodge::Wmemo = ();
  my $v = HodgePt::WP([], $J); my $r = scalar @$J; my $C = Comps::load("$FindBin::Bin/engine_h${h}_r${r}.txt");
  my $ok = 1; for my $k (0..$h) { my $a = Comps::comp_at($C, $h, $r, $k, $J); my $b = $v->[$k] // Math::BigRat->new(0); $ok = 0 if $a != $b; }
  print "h=$h J=(@$J): ", ($ok ? "agree" : "DISAGREE"), "  [k=0: ", $v->[0], "]\n"; }
