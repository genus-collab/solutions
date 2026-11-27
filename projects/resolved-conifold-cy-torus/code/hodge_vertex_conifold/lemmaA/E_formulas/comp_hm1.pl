# comp_hm1.pl h : the component next to lambda_g.  By Mumford's relation (Chern roots, prod_i (1 - alpha_i^2) = 1):
#   c(E^vee) e_{h-1}(alpha(1+alpha)) = lambda_{h-1} - h lambda_h + lambda_h Sum_{n odd} n! ch_n,
# so Xi^(h-1)_h(R) = X_a(R) - h Xi^(h)_h(R) + X_b(R) with X_a = Int K lambda_{h-1} prod lam, X_b = Int K lambda_h Sum n! ch_n prod lam.
# Checks the identity against the engine component and tests where X_a and X_b vanish separately (r = 1, 2).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use Hodge; use HodgePt;
sub R { Math::BigRat->new(shift) }
my $h = shift // 4; my $D = 3*$h - 3;
Hodge::setup($h);
my @lam = Hodge::lambdas($h, $D);
my $Xa = { map { $_ => $lam[$h-1]{$_} } keys %{$lam[$h-1]} };
my $sumch = {}; my @CHL = (1, 3, 5, 7); for my $i (0..$Hodge::NCH-1) { my $l = $CHL[$i]; my @e = (0) x $Hodge::NCH; $e[$i] = 1; my $f = R(1); $f *= $_ for 2..$l; $sumch->{join(',', @e)} = [$f]; }
my $Xb = Hodge::chmul($lam[$h], $sumch, $D);
sub with_class { my ($cls, $pts) = @_; local $Hodge::OM = $cls; %Hodge::Tmemo = (); %HodgePt::memo = (); my @v; for my $p (@$pts) { my $x = HodgePt::WP([], $p); push @v, $x->[0] // R(0); } %Hodge::Tmemo = (); %HodgePt::memo = (); return @v; }
my @pts; for my $j (1..2*$h+1) { push @pts, [$j]; } for my $a (2..2*$h) { for my $b ($a..2*$h) { next if $a + $b > 2*$h + 1; push @pts, [$a, $b]; } }
my @va = with_class($Xa, \@pts); my @vb = with_class($Xb, \@pts);
my @full; { %Hodge::Tmemo = (); %HodgePt::memo = (); for my $p (@pts) { my $x = HodgePt::WP([], $p); push @full, [$x->[$h-1] // R(0), $x->[$h] // R(0)]; } }
for my $i (0..$#pts) { my $p = $pts[$i]; my $D = 0; $D += $_ - 1 for @$p; my $pred = $va[$i] - $h * $full[$i][1] + $vb[$i];
  printf "h=%d R={%s} D=%d: Xi^(h-1) %s identity; X_a = %s, X_b = %s%s\n", $h, join(',', @$p), $D, ($pred == $full[$i][0] ? "=" : "!="), $va[$i], $vb[$i], ($D <= 2*$h - 3 ? "   [Lemma A range]" : ""); }
