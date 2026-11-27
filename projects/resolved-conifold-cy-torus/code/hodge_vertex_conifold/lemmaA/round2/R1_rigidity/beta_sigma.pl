#!/usr/bin/perl
# beta_sigma.pl -- rewrite beta_{k,L}(tau) = theta^{[L odd]} f_{k,L}(sigma), theta = tau + 1/2, sigma = tau(1+tau).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use r1lib;
my $file = shift // "$FindBin::Bin/out/beta_poly_K14.txt";
open my $fh, '<', $file or die; my %B;
while (<$fh>) { next unless /^beta (\d+) (\d+) : (.*)$/; $B{"$1,$2"} = [ map { R($_) } split ' ', $3 ]; } close $fh;
sub pdiv { my ($p, $d) = @_; my @p = @$p; my @q; my $dd = $#$d;   # polynomial division, returns (q, r)
  while (@p - 1 >= $dd && grep { !$_->is_zero } @p) { my $c = $p[-1] / $d->[-1]; my $sh = $#p - $dd; $q[$sh] = $c; for my $i (0..$dd) { $p[$sh+$i] = $p[$sh+$i] - $c*$d->[$i]; } pop @p; }
  @q = map { $_ // R(0) } @q; return (\@q, \@p); }
for my $key (sort { my @a = split /,/, $a; my @b = split /,/, $b; $a[0] <=> $b[0] || $a[1] <=> $b[1] } keys %B) {
  my ($k, $L) = split /,/, $key; my @p = @{$B{$key}}; next unless grep { !$_->is_zero } @p;
  my $pref = '';
  if ($L % 2) { my ($q, $r) = pdiv(\@p, [R('1/2'), R(1)]); die "not divisible by theta: $key" if grep { !$_->is_zero } @$r; @p = @$q; $pref = 'theta*'; }
  # express in sigma: repeatedly divide by sigma = tau + tau^2
  my @f; my $cur = [@p];
  while (grep { !$_->is_zero } @$cur) { my ($q, $r) = pdiv($cur, [R(0), R(1), R(1)]); my @r = @$r; push @r, R(0) while @r < 2;
    die "not symmetric: $key" unless $r[1]->is_zero; push @f, $r[0]; $cur = $q; last unless @$q; }
  # @f = coefficients of sigma^0, sigma^1, ...
  my $kap = undef;
  printf "f_{%d,%d} = %s[ %s ]\n", $k, $L, $pref, join(', ', map { "$_" } @f);
}
