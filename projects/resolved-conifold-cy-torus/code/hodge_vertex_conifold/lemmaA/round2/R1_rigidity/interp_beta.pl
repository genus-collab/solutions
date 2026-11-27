#!/usr/bin/perl
# interp_beta.pl -- the closed-form relation coefficients beta_{k,L}(tau) (relation_k = sum_L beta_{k,L} y_L, y_L = Pi_{L+1} Xi(L+1))
# as exact polynomials in tau: modular evaluation (kpmod.pm) at tau = 1..NT for several primes, Newton interpolation mod p,
# CRT and rational reconstruction.  One extra evaluation point per prime is used as a check of the degree bound.
# Usage: perl interp_beta.pl K [nprimes]     output: lines "beta k L : c0 c1 c2 ..." (coefficients of tau^0, tau^1, ...)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use kpmod; use Math::BigInt;
my ($K, $np) = @ARGV; $np //= 6;
my $NT = $ENV{NTPTS} // (2*$K + 2);   # tau-degree of beta_{k,L} is <= k-L+1 <= K (observed); NTPTS = K+4 suffices with checks
my @primes; { my $c = 2147483647; while (@primes < $np) { my $isp = 1; for (my $d = 3; $d*$d <= $c; $d += 2) { if ($c % $d == 0) { $isp = 0; last; } } push @primes, $c if $isp && $c % 2; $c -= 2; } }
my %res;   # res{"k,L"}[pi] = arrayref of poly coefficients mod p
for my $pi (0..$#primes) { set_prime($primes[$pi]); my $P = $primes[$pi];
  my @vals;   # vals[t] = rows
  for my $t (1..$NT) { my $S = kp_setup(md($t), 1, $K); my $rows = kp_rows($S); $vals[$t] = { map { ($_->[0] => $_->[1]) } @$rows }; }
  for my $k (1..$K) { for my $L (1..$k) {
      my @x = (1..$NT-1); my @y = map { $vals[$_]{$k}[$L-1] } @x;
      # Newton divided differences mod p
      my @dd = @y; my $n = @x;
      for my $j (1..$n-1) { for (my $i = $n-1; $i >= $j; $i--) { $dd[$i] = mm(md($dd[$i] - $dd[$i-1]), mi(md($x[$i] - $x[$i-$j]))); } }
      # convert Newton form to monomial coefficients
      my @c = (0) x $n; my @basis = (1);   # basis poly prod (t - x_i)
      for my $j (0..$n-1) { for my $i (0..$#basis) { $c[$i] = ($c[$i] + mm($dd[$j], $basis[$i])) % $P; }
        my @nb = (0) x (@basis + 1); for my $i (0..$#basis) { $nb[$i+1] = ($nb[$i+1] + $basis[$i]) % $P; $nb[$i] = md($nb[$i] - mm($x[$j], $basis[$i])); } @basis = @nb; }
      # check at the last point
      my $tc = $NT; my $v = 0; my $pw = 1; for my $i (0..$n-1) { $v = ($v + mm($c[$i], $pw)) % $P; $pw = mm($pw, $tc); }
      die "degree check failed k=$k L=$L p=$P" unless $v == $vals[$tc]{$k}[$L-1];
      $res{"$k,$L"}[$pi] = \@c; } } }
# CRT + rational reconstruction
my $M = Math::BigInt->new(1); $M *= $_ for @primes;
sub crt { my @r = @_; my $x = Math::BigInt->new(0);
  for my $i (0..$#primes) { my $p = $primes[$i]; my $Mi = $M->copy->bdiv($p); my $inv = $Mi->copy->bmodinv($p); $x += Math::BigInt->new($r[$i]) * $Mi * $inv; } return $x->bmod($M); }
sub ratrec { my $r = shift; my ($a0, $a1) = ($M->copy, $r->copy); my ($b0, $b1) = (Math::BigInt->new(0), Math::BigInt->new(1)); my $bound = $M->copy->bdiv(2)->bsqrt;
  while ($a1 > $bound) { my $q = $a0->copy->bdiv($a1); ($a0, $a1) = ($a1, $a0 - $q*$a1); ($b0, $b1) = ($b1, $b0 - $q*$b1); }
  die "ratrec failed" if $b1->copy->babs > $bound; my ($n, $d) = ($a1, $b1); if ($d < 0) { $n = -$n; $d = -$d; } my $g = Math::BigInt::bgcd($n, $d); return ($n/$g) . "/" . ($d/$g); }
for my $k (1..$K) { for my $L (1..$k) { my $nc = scalar @{$res{"$k,$L"}[0]}; my @out; my $last = -1;
    for my $i (0..$nc-1) { my $r = crt(map { $res{"$k,$L"}[$_][$i] } 0..$#primes); my $q = ratrec($r); $q =~ s{/1$}{}; push @out, $q; $last = $i if $q ne '0'; }
    print "beta $k $L : ", join(' ', @out[0..($last < 0 ? 0 : $last)]), "\n"; } }
