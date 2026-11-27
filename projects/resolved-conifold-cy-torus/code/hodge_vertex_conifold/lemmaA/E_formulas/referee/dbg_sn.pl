# dbg_sn.pl : debugging aid for t_sn.pl (compare slice data of the two engines at genus 2, m = 4)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use lib "$FindBin::Bin/.."; use Math::BigRat; use Hodge; use HodgePt;
sub Q { Math::BigRat->new(shift) }
Hodge::setup(2);
my @C; open my $fh, '<', "$FindBin::Bin/ref_engine_h2_r2.txt" or die; while (my $line = <$fh>) { chomp $line; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4;
  for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($mm, $cf) = split /:/, $t; push @{$C[$k]}, [[split /,/, $mm], Q($cf)]; } } close $fh;
for my $j (0..6) { my $v = HodgePt::WP([], [4, $j]); my @mine;
  for my $k (0..2) { my $s = Q(0); for my $t (@{$C[$k]}) { my ($e, $cf) = @$t; $s += $cf * Q(4) ** $e->[0] * Q($j) ** $e->[1]; } push @mine, $s; }
  print "j=$j author: ", join(' ', map { $_ // 0 } @$v), " | referee: @mine\n"; }
