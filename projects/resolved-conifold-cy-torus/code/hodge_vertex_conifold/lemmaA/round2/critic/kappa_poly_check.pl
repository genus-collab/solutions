#!/usr/bin/perl
# Critic check of brief s.2(f): kappa_D(tau) = Sum_{T at level D, parts >= 2} (1/|Aut T|) (J-1)...(J-r+1) prod_t (-1)^t Pi_t(tau)/t
# equals [y^D](tau y + log((e^y-1)/y)), i.e. tau + 1/2 at D = 1 and B_D/(D D!) for D >= 2.
# kappa_D is a polynomial in tau of degree <= D (deg Pi_t = t-1), so agreement at D+1 distinct tau proves the identity at that D.
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
sub Pi { my ($t,$tau)=@_; my $p=R(1); for my $a (1..$t-1){ $p*= ($tau*$t + $a); $p/= $a; } $p }
sub fact { my $n=shift; my $f=R(1); $f*=$_ for 1..$n; $f }
sub parts { my ($n,$max)=@_; return ([]) if $n==0; my @o; for my $k (reverse 1..($max<$n?$max:$n)){ push @o, [$k,@$_] for parts($n-$k,$k); } @o }
my $DMAX = shift // 14;
my @B=(R(1)); for my $m (1..$DMAX+1){ my $s=R(0); for my $k (0..$m-1){ $s+= fact($m+1)/(fact($k)*fact($m+1-$k))*$B[$k]; } $B[$m]= -$s/($m+1); }
my $bad=0; my $n=0;
for my $D (1..$DMAX) {
  my @lams = parts($D,$D);
  my @taus = map { R($_) } (0..$D);                       # D+1 integer points
  push @taus, R("3/11"), R("-7/5");                       # two extra non-integer points
  for my $tau (@taus) {
    my $sum=R(0);
    for my $lam (@lams) {
      my @t = map { $_+1 } @$lam; my $r=@t; my $J=$D+$r;
      my %m; $m{$_}++ for @t; my $aut=R(1); $aut*=fact($_) for values %m;
      my $w = fact($J-1)/fact($J-$r)/$aut;
      for my $tt (@t){ $w *= ( ($tt%2) ? -1 : 1 ) * Pi($tt,$tau) / $tt; }
      $sum += $w;
    }
    my $pred = $D==1 ? $tau + R("1/2") : $B[$D]/($D*fact($D));
    $n++; if ($sum != $pred) { $bad++; print "MISMATCH D=$D tau=$tau kappa=$sum pred=$pred\n"; }
  }
  printf "D=%2d: %3d multisets, %2d tau values, kappa = %s\n", $D, scalar(@lams), scalar(@taus), ($D==1 ? "tau+1/2" : "".($B[$D]/($D*fact($D))));
}
print "total $n evaluations, $bad mismatches\n";
