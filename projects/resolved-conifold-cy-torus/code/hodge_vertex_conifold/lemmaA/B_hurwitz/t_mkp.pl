#!/usr/bin/perl
# Exploration (outlook, CY extension): test of the first "modified KP" equation linking the charge-0 and charge-1 tau functions
# of the 2D Toda sequence of Okounkov (2000), which involves only derivatives in the first set of times p:
#   s * 2 d/dp_2 (F - F1) + d^2/dp_1^2 (F + F1) + (d/dp_1 (F - F1))^2 = 0,   F = H(p), F1 = H(e^b o p),  s = +1 or -1,
# H = connected single Hurwitz series as in t_toda.pl.  (t_2 = p_2/2.)
use strict; use warnings; use lib '.'; use HB; use HBconn; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my $DMAX = shift // 6; my $BMAX = shift // 8;
my %H;
for my $d (1..$DMAX) { for my $mu (@{HB::partitions($d)}) { my @nu = grep { $_ >= 2 } @$mu; my $key = join(",", @$mu);
  my %aut; $aut{$_}++ for @nu; my $autnu = Math::BigInt->new(1); $autnu *= HB::fact($_) for values %aut;
  my $C = HBconn::conn_counts(\@nu, $d, $BMAX + 2);
  for my $b (0..$BMAX + 2) { next unless exists $C->{$d}{$b}; my $v = BR($C->{$d}{$b}) / (BR(HB::fact($d)) * BR($autnu) * BR(HB::fact($b)));
    $H{$key}{$b} = $v unless $v == 0 } } }
sub deg { my $k = shift; my $s = 0; $s += $_ for split /,/, $k; return $s }
sub mono_mul { my ($ma, $mb) = @_; my @x = ((split /,/, $ma), (split /,/, $mb)); @x = sort { $b <=> $a } @x; return join(",", @x) }
sub smul { my ($A, $B, $dmax, $bmax) = @_; my %P;
  for my $ka (keys %$A) { my $da = deg($ka); for my $kb (keys %$B) { next if $da + deg($kb) > $dmax; my $k = mono_mul($ka, $kb);
    for my $ba (keys %{$A->{$ka}}) { for my $bb (keys %{$B->{$kb}}) { next if $ba + $bb > $bmax;
      $P{$k}{$ba + $bb} = ($P{$k}{$ba + $bb} // BR(0)) + $A->{$ka}{$ba} * $B->{$kb}{$bb} } } } }
  return \%P }
sub sadd { my ($A, $B, $c) = @_; my %P; for my $k (keys %$A) { for my $b (keys %{$A->{$k}}) { $P{$k}{$b} = $A->{$k}{$b}->copy } }
  for my $k (keys %$B) { for my $b (keys %{$B->{$k}}) { $P{$k}{$b} = ($P{$k}{$b} // BR(0)) + $c * $B->{$k}{$b} } } return \%P }
# F1 = H(e^b o p): multiply p_mu b^B by e^{|mu| b}
my %F1;
for my $k (keys %H) { my $d = deg($k); for my $b (keys %{$H{$k}}) { for (my $e = 0; $b + $e <= $BMAX + 2; $e++) {
  $F1{$k}{$b + $e} = ($F1{$k}{$b + $e} // BR(0)) + $H{$k}{$b} * BR(Math::BigInt->new($d) ** $e) / BR(HB::fact($e)) } } }
sub dpk { my ($A, $kk) = @_; my %P; # d/dp_kk
  for my $key (keys %$A) { my @mu = split /,/, $key; my $m = grep { $_ == $kk } @mu; next unless $m;
    my @rest = @mu; for my $i (0..$#rest) { if ($rest[$i] == $kk) { splice(@rest, $i, 1); last } } my $kr = join(",", @rest);
    for my $b (keys %{$A->{$key}}) { $P{$kr}{$b} = ($P{$kr}{$b} // BR(0)) + $A->{$key}{$b} * $m } }
  return \%P }
my $Dm = sadd(\%H, \%F1, BR(-1)); my $Sp = sadd(\%H, \%F1, BR(1));
my $T2 = dpk($Dm, 2); my $T11 = dpk(dpk($Sp, 1), 1); my $T1 = dpk($Dm, 1); my $Q = smul($T1, $T1, $DMAX - 2, $BMAX);
for my $s (1, -1) { my $tot = sadd(sadd($T11, $Q, BR(1)), $T2, BR(2 * $s)); my ($n, $nz) = (0, 0);
  for my $k (keys %$tot) { next if deg($k) > $DMAX - 2; for my $b (keys %{$tot->{$k}}) { next if $b > $BMAX; $n++; $nz++ if $tot->{$k}{$b} != 0 } }
  print "sign $s: $n coefficients (deg <= ", $DMAX - 2, ", beta^<= $BMAX), nonzero: $nz\n" }
