#!/usr/bin/perl
# Exact check of BS v2 r_l: [exp(z d1 d2/3) exp(s1^3/z) g(s)]_{s1=s2=0} = -sqrt(3) i/s3^2 * sum_l r_l (z/s3)^l
# with g = 9/(s3 sqrt(Disc)) (v2 measure) and, as control, g_old = sqrt(Disc)/s3^3 (old measure). s3 = 1.
use strict; use warnings; use Math::BigRat;
my $DEG = 12; my $LMAX = 4;
sub R { Math::BigRat->new(@_) }
sub mul { my ($p,$q)=@_; my %r; for my $a (keys %$p){ my ($i,$j)=split /,/,$a; for my $b (keys %$q){ my ($k,$l)=split /,/,$b; next if $i+$k>$DEG||$j+$l>$DEG; my $key=($i+$k).",".($j+$l); $r{$key}=(exists $r{$key}?$r{$key}:R(0)) + $p->{$a}*$q->{$b}; } } return \%r; }
sub add { my ($p,$q,$c)=@_; my %r=%$p; for my $b (keys %$q){ $r{$b}=(exists $r{$b}?$r{$b}:R(0)) + $c*$q->{$b}; } return \%r; }
# u = (s1^2 s2^2 - 4 s1^3 + 18 s1 s2 - 4 s2^3)/27  (Disc = -27 (1 - u) at s3 = 1)
my $u = { "2,2"=>R(1,27), "3,0"=>R(-4,27), "1,1"=>R(18,27), "0,3"=>R(-4,27) };
for my $mode ("v2","old") {
  # v2: g = 9/(sqrt(-27) sqrt(1-u)) = (-sqrt3 i) (1-u)^(-1/2);  old: g = sqrt(-27) sqrt(1-u) = (3 sqrt3 i)(1-u)^(1/2) = (-sqrt3 i)*(-3)(1-u)^(1/2)
  my $ex = $mode eq "v2" ? R(-1,2) : R(1,2); my $pref = $mode eq "v2" ? R(1) : R(-3);
  my %G = ("0,0"=>R(1)); my $up = {"0,0"=>R(1)}; my $bin = R(1);
  for my $k (1..$DEG) { $bin = $bin*($ex - Math::BigRat->new($k-1))/R($k); $up = mul($up,$u); %G = %{ add(\%G,$up,$bin*R((-1)**$k)) }; }
  # sum over n, m: (z/3)^n n! [s1^n s2^n] (s1^{3m}/(m! z^m) G) -> z^{n-m}
  my @r = map { R(0) } 0..$LMAX;
  for my $n (0..$DEG) { for my $m (0..int($n/3)) { my $l=$n-$m; next if $l>$LMAX; my $key=($n-3*$m).",".$n; next unless exists $G{$key};
      my $f=R(1); $f*=R($_) for 1..$n; my $mf=R(1); $mf*=R($_) for 1..$m;
      $r[$l] += $pref * $G{$key} * $f / ($mf * R(3)**$n); } }
  print "$mode measure: ", join("  ", map {"r_$_ = $r[$_]"} 0..$LMAX), "\n";
}
