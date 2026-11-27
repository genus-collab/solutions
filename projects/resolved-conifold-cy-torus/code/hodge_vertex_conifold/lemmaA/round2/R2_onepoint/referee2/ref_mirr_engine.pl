#!/usr/bin/perl
# ref_mirr_engine.pl TABLE h r : TR-free test of Lemma A (all components) from an engine table: every component Xi^(k)_h(R) must vanish
# for every multiset R of r positive parts with D(R) = |R| - r <= 2h-3 (parts >= 1, which includes M_irr(h) for parts >= 2).
# Also lists the first non-vanishing layer D = 2h-2 (should be non-zero somewhere).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($file, $h, $r) = @ARGV; my $C = RefR2::load_table($file, $h, $r);
my @tuples; my $rec; $rec = sub { my ($pre, $max, $left) = @_; if (@$pre == $r) { push @tuples, [@$pre]; return; } for my $x (reverse 1..($max < $left ? $max : $left)) { $rec->([@$pre, $x], $x, $left - $x); } };
$rec->([], 2*$h - 2 + $r, 2*$h - 2 + $r);
my ($n, $bad, $layer, $layernz) = (0, 0, 0, 0);
for my $t (@tuples) { my $J = 0; $J += $_ for @$t; my $D = $J - $r;
  for my $k (sort { $a <=> $b } keys %$C) { my $s = Q(0);
    for my $key (keys %{$C->{$k}}) { my @a = split /,/, $key; my $x = $C->{$k}{$key}->copy; $x *= Q($t->[$_])**$a[$_] for 0..$#a; $s += $x; }
    if ($D <= 2*$h - 3) { $n++; if (!$s->is_zero) { $bad++; print "NONZERO h=$h R=(@$t) k=$k: $s\n"; } }
    elsif ($D == 2*$h - 2) { $layer++; $layernz++ unless $s->is_zero; } } }
printf "genus %d r=%d: %d (R, component) values in the range D <= %d, %d non-zero; first layer D = %d: %d of %d component values non-zero\n",
  $h, $r, $n, 2*$h - 3, $bad, 2*$h - 2, $layernz, $layer;
