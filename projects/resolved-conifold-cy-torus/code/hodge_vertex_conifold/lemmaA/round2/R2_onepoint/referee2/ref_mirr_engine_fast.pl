#!/usr/bin/perl
# ref_mirr_engine_fast.pl TABLE h r : same test as ref_mirr_engine.pl (every component Xi^(k)_h(R) vanishes for D(R) <= 2h-3, parts >= 1;
# first layer D = 2h-2 reported), with the powers j_i^a precomputed once per multiset R (faster for large r = 3 tables).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefR2; use Math::BigRat;
sub Q { Math::BigRat->new(shift) }
my ($file, $h, $r) = @ARGV; my $C = RefR2::load_table($file, $h, $r);
# pre-split the table
my %T; my $amax = 0;
for my $k (keys %$C) { for my $key (keys %{$C->{$k}}) { my @a = split /,/, $key; for (@a) { $amax = $_ if $_ > $amax } push @{$T{$k}}, [\@a, $C->{$k}{$key}]; } }
my @tuples; my $rec; $rec = sub { my ($pre, $max, $left) = @_; if (@$pre == $r) { push @tuples, [@$pre]; return; } for my $x (reverse 1..($max < $left ? $max : $left)) { $rec->([@$pre, $x], $x, $left - $x); } };
$rec->([], 2*$h - 2 + $r, 2*$h - 2 + $r);
my ($n, $bad, $layer, $layernz) = (0, 0, 0, 0);
for my $t (@tuples) { my $J = 0; $J += $_ for @$t; my $D = $J - $r; next if $D > 2*$h - 2;
  my @pw; for my $i (0..$r-1) { my @p = (Q(1)); my $x = Q($t->[$i]); push @p, $p[-1] * $x for 1..$amax; push @pw, \@p; }
  for my $k (sort { $a <=> $b } keys %T) { my $s = Q(0);
    for my $term (@{$T{$k}}) { my ($a, $c) = @$term; my $x = $c->copy; $x *= $pw[$_][$a->[$_]] for 0..$r-1; $s += $x; }
    if ($D <= 2*$h - 3) { $n++; if (!$s->is_zero) { $bad++; print "NONZERO h=$h R=(@$t) k=$k: $s\n"; } }
    else { $layer++; $layernz++ unless $s->is_zero; } } }
printf "genus %d r=%d: %d (R, component) values in the range D <= %d, %d non-zero; first layer D = %d: %d of %d component values non-zero\n",
  $h, $r, $n, 2*$h - 3, $bad, 2*$h - 2, $layernz, $layer;
