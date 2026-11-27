#!/usr/bin/perl
# Numerical test of the Toda equation for single Hurwitz numbers (Okounkov 2000, specialised to p' = (1,0,0,...)):
#   d/dp_1 ( E H ) = exp( H(e^b o p) - 2 H(p) + H(e^-b o p) ),   E = sum_k k p_k d/dp_k,  (e^b o p)_k = e^(k b) p_k,
# where H(p; b) = sum_{g,mu} H_{g,mu} p_mu b^B / B!  is the connected Hurwitz series (B = number of simple branch points).
# Truncation: p-degree (= covering degree) <= DMAX on the left (<= DMAX-1 after d/dp_1), beta-order <= BMAX.
use strict; use warnings; use lib '.'; use HB; use HBconn; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my $DMAX = shift // 5; my $BMAX = shift // 7;
# series: { "mu" => { b => coef } } with mu a sorted (descending) comma list ("" for the empty partition); coef of p_mu beta^b
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
# Delta H: multiply p_mu beta^b by 2cosh(|mu| beta) - 2 = sum_{k>=1} 2 |mu|^(2k) beta^(2k)/(2k)!
my %D;
for my $k (keys %H) { my $d = deg($k); for my $b (keys %{$H{$k}}) { for (my $e = 2; $b + $e <= $BMAX; $e += 2) {
  $D{$k}{$b + $e} = ($D{$k}{$b + $e} // BR(0)) + $H{$k}{$b} * 2 * BR(Math::BigInt->new($d) ** $e) / BR(HB::fact($e)) } } }
# exp(D) truncated at degree DMAX-1
my %EXP = ("" => { 0 => BR(1) }); my $term = { "" => { 0 => BR(1) } };
for my $n (1..$DMAX) { $term = smul($term, \%D, $DMAX - 1, $BMAX); my %t2; for my $k (keys %$term) { for my $b (keys %{$term->{$k}}) { $t2{$k}{$b} = $term->{$k}{$b} / $n } } $term = \%t2;
  %EXP = %{ sadd(\%EXP, $term, BR(1)) } }
# LHS: d/dp_1 E H
my %L;
for my $k (keys %H) { my @mu = split /,/, $k; my $m1 = grep { $_ == 1 } @mu; next unless $m1; my $d = deg($k);
  my @rest = @mu; for my $i (0..$#rest) { if ($rest[$i] == 1) { splice(@rest, $i, 1); last } } my $kr = join(",", @rest);
  for my $b (keys %{$H{$k}}) { next if $b > $BMAX; $L{$kr}{$b} = ($L{$kr}{$b} // BR(0)) + $H{$k}{$b} * $d * $m1 } }
my ($n, $bad) = (0, 0);
my %keys; $keys{"$_"} = 1 for keys %L; $keys{"$_"} = 1 for keys %EXP;
for my $k (sort keys %keys) { next if deg($k) > $DMAX - 1; for my $b (0..$BMAX) { my $l = $L{$k}{$b} // BR(0); my $r = $EXP{$k}{$b} // BR(0);
  next if $l == 0 && $r == 0; $n++; if ($l != $r) { $bad++; print "MISMATCH p_[$k] beta^$b: LHS $l RHS $r\n" } } }
print "Toda test (deg <= $DMAX, beta^<= $BMAX): $n nonzero coefficients compared, $bad mismatches\n";
