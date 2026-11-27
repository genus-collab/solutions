# t_e.pl h1 h2 ... : the author's conjectural closed formula 0(e) for the component next to lambda_g (two leaves), tested on the
# REFEREE polynomials ref_engine_h<h>_r2.txt (genus 7 was never seen by the author: a parameter-free prediction):
#   (-1)^{h-1} Xi^(h-1)_h(j1,j2) = prod_{l=1}^{2h-1}(J-l) [l^{2h}] (l S'(l) - J S(l) log S(l))
#        + Sum_{k=1}^{h} (j1)_k (j2)_k prod_{l=2k}^{2h-1}(J-l) [l^{2h}] S(l) G_k(l),
#   G_k = (-1)^{k+1} (1/k) Sum_n C(2n-k, k-1) c_n l^{2n},  log S = Sum_n c_n l^{2n},  S(l) = (l/2)/sinh(l/2).
# Implemented independently (log S from the sinh series; c_n compared with -B_{2n}/(2n (2n)!)).  Identity check on the full grid
# 0 <= j1, j2 <= 3h+1 (a polynomial of degree <= 3h-1 vanishing on it is zero).
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use RefHodge;
sub Q { Math::BigRat->new(shift) }
sub fact { my $n = shift; my $x = Q(1); $x *= $_ for 2..$n; $x }
sub binom { my ($n, $k) = @_; return Q(0) if $k < 0 || $k > $n || $n < 0; my $x = Q(1); for my $i (1..$k) { $x = $x * ($n - $k + $i) / $i } $x }
for my $h (@ARGV) { my $K = $h;
  my @sh = map { Q(1) / (Q(4) ** $_ * fact(2*$_ + 1)) } 0..$K;
  my @lg = (Q(0)); for my $n (1..$K) { my $s = $sh[$n] * $n; for my $m (1..$n-1) { $s -= $m * $lg[$m] * $sh[$n - $m]; } push @lg, $s / $n; }
  my @c = map { -$_ } @lg;                                   # log S
  for my $n (1..$K) { die "c_n mismatch" unless $c[$n] == -RefHodge::bern(2*$n) / (2*$n * fact(2*$n)); }
  my @S = (Q(1)); for my $n (1..$K) { my $s = Q(0); for my $m (1..$n) { $s += $m * $c[$m] * $S[$n - $m]; } $S[$n] = $s / $n; }
  my $smul = sub { my ($a, $b) = @_; my @r = map { Q(0) } 0..$K; for my $i (0..$K) { for my $j (0..$K-$i) { $r[$i+$j] += $a->[$i] * $b->[$j]; } } \@r };
  my $SlogS = $smul->(\@S, \@c);
  my @co; for my $k (1..$h) { my @g = map { Q(0) } 0..$K; for my $n ($k..$K) { $g[$n] = ((($k + 1) % 2) ? Q(-1) : Q(1)) / $k * binom(2*$n - $k, $k - 1) * $c[$n]; } $co[$k] = $smul->(\@S, \@g)->[$h]; }
  # data
  my %P; open my $fh, '<', "$FindBin::Bin/ref_engine_h${h}_r2.txt" or die; while (my $line = <$fh>) { chomp $line; my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4; next unless $k == $h - 1;
    for my $t (split /\s*;\s*/, $rest // '') { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; $P{$m} = Q($cf); } } close $fh;
  my $sgn = (($h - 1) % 2) ? Q(-1) : Q(1);
  my ($n, $bad) = (0, 0);
  for my $j1 (0..3*$h + 1) { for my $j2 ($j1..3*$h + 1) { my $J = $j1 + $j2;
      my $data = Q(0); for my $m (keys %P) { my ($a, $b) = split /,/, $m; $data += $P{$m} * Q($j1) ** $a * Q($j2) ** $b; } $data *= $sgn;
      my $pr = Q(1); $pr *= ($J - $_) for 1..(2*$h - 1); my $pred = $pr * (2*$h * $S[$h] - $J * $SlogS->[$h]);
      for my $k (1..$h) { my $ff = Q(1); for my $i (0..$k-1) { $ff *= ($j1 - $i) * ($j2 - $i); } next if $ff == 0; my $t = $ff * $co[$k]; $t *= ($J - $_) for (2*$k)..(2*$h - 1); $pred += $t; }
      $n++; $bad++ if $pred != $data; } }
  printf "genus %d: closed formula 0(e) for (-1)^{h-1} Xi^(h-1)(j1,j2) vs referee engine on the grid 0..%d: %d points, %d mismatches\n", $h, 3*$h + 1, $n, $bad; }
