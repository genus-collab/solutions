#!/usr/bin/perl
# kappa_D(tau) = Sum_{T at level D} (1/|Aut T|) (J-1)...(J-r+1) prod_i (-1)^{t_i} Pi_{t_i}(tau)/t_i
# i.e. Route D's relation T_D after substituting Route A's Corollary E (all level-D values in terms of Xi(D+1)),
# divided by Pi_{D+1} Xi(D+1). Prediction: kappa_D = [y^D] (tau y + log((e^y - 1)/y)).
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(@_) }
sub Pi { my ($t,$tau)=@_; my $p=R(1); for my $a (1..$t-1){ $p*= ($tau*$t + $a); } for my $a (1..$t-1){ $p/= $a; } $p }
sub fact { my $n=shift; my $f=R(1); $f*=$_ for 1..$n; $f }
sub parts { my ($n,$max)=@_; return ([]) if $n==0; my @o; for my $k (reverse 1..($max<$n?$max:$n)){ push @o, [$k,@$_] for parts($n-$k,$k); } @o }
# Bernoulli numbers
my @B=(R(1)); for my $m (1..20){ my $s=R(0); for my $k (0..$m-1){ my $c=fact($m+1)/(fact($k)*fact($m+1-$k)); $s+= $c*$B[$k]; } $B[$m]= -$s/($m+1); }
for my $tau (R(0), R("2/7"), R("3/2"), R("-5/3")) {
  print "tau = $tau\n";
  for my $D (1..10) {
    my $sum=R(0);
    for my $lam (parts($D,$D)) {
      my @t = map { $_+1 } @$lam; my $r=@t; my $J=$D+$r;
      my %m; $m{$_}++ for @t; my $aut=R(1); $aut*=fact($_) for values %m;
      my $w = fact($J-1)/fact($J-$r)/$aut;
      for my $t (@t){ $w *= ( ($t%2) ? -1 : 1 ) * Pi($t,$tau) / $t; }
      $sum += $w;
    }
    my $pred = $D==1 ? $tau + R("1/2") : ($D%2 ? R(0) : $B[$D]/($D*fact($D)));
    printf "  D=%2d  kappa=%-28s pred=%-28s %s\n", $D, "$sum", "$pred", ($sum==$pred ? "OK" : "MISMATCH");
  }
}
