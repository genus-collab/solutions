package ratrec;
# CRT over several primes and rational reconstruction (Math::BigInt), plus a list of primes below 2^31.
use strict; use warnings; use Math::BigInt; use Math::BigRat;
use Exporter 'import';
our @EXPORT_OK = qw(primes_below crt ratrecon recon);
sub isprime { my $n=shift; return 0 if $n<2; return 1 if $n<4; return 0 if $n%2==0; for (my $d=3; $d*$d<=$n; $d+=2){ return 0 if $n%$d==0; } 1 }
sub primes_below { my ($B,$k)=@_; my @p; my $n=$B-1; while (@p<$k){ push @p,$n if isprime($n); $n--; } @p }
sub crt { my ($res,$mods)=@_; my $M=Math::BigInt->new(1); my $x=Math::BigInt->new(0);
  for my $i (0..$#$res){ my $m=Math::BigInt->new($mods->[$i]); my $r=Math::BigInt->new($res->[$i]);
    # x = x + M * ((r - x) * M^-1 mod m)
    my $t = ($r - $x) % $m; my $Mi = $M->copy->bmodinv($m); $t = ($t * $Mi) % $m; $x = $x + $M*$t; $M = $M*$m; }
  ($x % $M, $M) }
sub ratrecon { my ($a,$M)=@_; # find n/d with |n|,d <= sqrt(M/2), n = a d mod M
  my $B = ($M/2)->bsqrt; my ($r0,$r1)=($M->copy,$a->copy % $M); my ($t0,$t1)=(Math::BigInt->new(0),Math::BigInt->new(1));
  while ($r1 > $B){ my $q=$r0/$r1; ($r0,$r1)=($r1,$r0-$q*$r1); ($t0,$t1)=($t1,$t0-$q*$t1); }
  return undef if $t1->copy->babs > $B || $t1->is_zero;
  my $n=$r1; my $d=$t1; if ($d<0){ $n=-$n; $d=-$d; } Math::BigRat->new("$n/$d") }
sub recon { my ($res,$mods)=@_; my ($x,$M)=crt($res,$mods); ratrecon($x,$M) }
1;
