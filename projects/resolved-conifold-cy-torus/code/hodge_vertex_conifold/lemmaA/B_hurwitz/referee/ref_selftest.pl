#!/usr/bin/perl
# Referee self-tests of RefHur.pm: (1) cut-and-join versus characters for Zhat; (2) brute force in S_4/S_5;
# (3) genus 0 (Hurwitz formula) and genus 1; (4) pure-unit connected series versus DYZ's closed forms H_2, H_3
# (Moscow Math. J. 17 (2017), eq. (33));  (5) Lambert structure and Prop D agree.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefHur; use Math::BigRat; use Math::BigInt;
sub BR { Math::BigRat->new(shift) } sub BI { Math::BigInt->new(shift) }
my ($n, $bad) = (0, 0); sub rep { my ($w, $a, $b) = @_; $n++; my $ok = ($a == $b); $bad++ unless $ok; print "$w: $a vs $b ", ($ok ? "ok" : "MISMATCH"), "\n" }
# (1) cut-and-join vs characters
for my $d (1..9) { for my $mu (@{RefHur::parts_of($d)}) { my @nu = grep { $_ >= 2 } @$mu; my $c = grep { $_ == 1 } @$mu;
  my $a = RefHur::zhat_cj(\@nu, $c, 12); my $b = RefHur::zhat_char(\@nu, $c, 12); my $same = 1; for (0..12) { $same = 0 unless $a->[$_] == $b->[$_] }
  $n++; $bad++ unless $same; print "Zhat mu=(@$mu) b<=12: cut-and-join vs characters ", ($same ? "ok" : "MISMATCH"), "\n" if !$same || $d <= 4 } }
# (2) brute force: number of b-tuples of transpositions in S_n with product equal to a fixed pi
sub pmul { my ($p, $q) = @_; [ map { $p->[$q->[$_]] } 0..$#$q ] }
sub brute { my ($n, $b, $pi) = @_; my @tr; for my $i (0..$n - 1) { for my $j ($i + 1..$n - 1) { my @p = 0..$n - 1; @p[$i, $j] = @p[$j, $i]; push @tr, \@p } }
  my %cnt = (join(",", 0..$n - 1) => 1);
  for (1..$b) { my %new; for my $k (keys %cnt) { my $p = [split /,/, $k]; for my $t (@tr) { $new{join(",", @{pmul($p, $t)})} += $cnt{$k} } } %cnt = %new }
  return $cnt{join(",", @$pi)} // 0 }
for my $case ([4, 3, [1,0,2,3]], [4, 4, [1,2,0,3]], [4, 4, [1,0,3,2]], [4, 5, [1,2,3,0]], [5, 4, [1,2,0,3,4]], [5, 4, [1,2,3,4,0]], [5, 5, [1,0,2,3,4]], [5, 6, [1,2,3,4,0]]) {
  my ($nn, $bb, $pi) = @$case; my %seen; my @type; for my $s (0..$nn - 1) { next if $seen{$s}; my $l = 0; my $x = $s; while (!$seen{$x}) { $seen{$x} = 1; $x = $pi->[$x]; $l++ } push @type, $l }
  @type = reverse sort { $a <=> $b } @type; my $T = RefHur::cj_table($nn, $bb)->{join(",", @type)}[$bb];
  rep("N_$bb(type @type) in S_$nn: cut-and-join vs brute force", $T, BI(brute($nn, $bb, $pi))) }
# (3) genus 0 and genus 1 values
for my $R ([2], [3], [2,3], [2,2,4]) { my $J = 0; $J += $_ for @$R; my $r = @$R; my $A = RefHur::A_list($R, 0, 5, 'char');
  for my $c (0..5) { next if $r + $c < 3; rep("A_0(@$R,1^$c) vs (J+c)^(r+c-3)", $A->[$c], BR(BI($J + $c) ** ($r + $c - 3))) } }
for my $j (2..6) { my $A = RefHur::A_list([$j], 1, 1, 'char'); rep("A_1($j) vs (j-1)/24", $A->[0], BR($j - 1) / 24); rep("A_1($j,1) vs j^2/24", $A->[1], BR($j * $j) / 24) }
# (4) pure units: G_{g,empty}(z) = sum_c A_g(1^c) z^c/c!  versus DYZ (33): H_2 = (T^2+6T^3)/(1440(1-T)^5), H_3 = (9T^2+548T^3+3482T^4+3816T^5+720T^6)/(725760(1-T)^10), z = T e^{-T}
sub Tser { my $N = shift; [ BR(0), map { BR(BI($_) ** ($_ - 1)) / BR(RefHur::fac($_)) } 1..$N ] }   # T(z) = sum d^{d-1} z^d/d!
sub smul { my ($a, $b, $N) = @_; my @c = map { BR(0) } 0..$N; for my $i (0..$N) { for my $j (0..$N - $i) { $c[$i + $j] += $a->[$i] * $b->[$j] } } \@c }
sub spow { my ($a, $k, $N) = @_; my $r = [BR(1), map { BR(0) } 1..$N]; $r = smul($r, $a, $N) for 1..$k; $r }
my $N = 8; my $T = Tser($N); my $omT = [ map { -$_ } @$T ]; $omT->[0] += 1;
my $inv = [BR(1), map { BR(0) } 1..$N]; { my @c = (BR(1)); for my $m (1..$N) { my $s = BR(0); for my $k (1..$m) { $s += $omT->[$k] * $c[$m - $k] } push @c, -$s } $inv = \@c }
my %dyz = (2 => [ [2, 1], [3, 6] ], 3 => [ [2, 9], [3, 548], [4, 3482], [5, 3816], [6, 720] ]); my %den = (2 => [1440, 5], 3 => [725760, 10]);
for my $g (2, 3) { my $num = [ map { BR(0) } 0..$N ]; for my $tc (@{$dyz{$g}}) { my $p = spow($T, $tc->[0], $N); $num->[$_] += BR($tc->[1]) * $p->[$_] for 0..$N }
  my $ser = smul($num, spow($inv, $den{$g}[1], $N), $N); $ser = [ map { $_ / BR($den{$g}[0]) } @$ser ];
  my $A = RefHur::A_list([], $g, $N, 'char');
  for my $c (1..$N) { rep("pure-unit genus $g, z^$c: Hurwitz (referee engine) vs DYZ closed form", $A->[$c] / BR(RefHur::fac($c)), $ser->[$c]) } }
# (5) Lambert structure (overdetermined) and Prop D, small cases; R = empty gives chi(M_g) = B_2g/(2g(2g-2))
for my $case ([[], 2], [[], 3], [[2], 1], [[3], 2], [[2,2], 2], [[2,3], 1], [[4], 3]) { my ($R, $h) = @$case;
  my ($p0, $ok, $pd) = RefHur::xi0($R, $h, 3, 'char'); $n++; $bad++ unless $ok; print "Lambert structure h=$h R=(@$R): ", ($ok ? "ok" : "FAILED"), "\n";
  rep("p0 vs Prop D, h=$h R=(@$R)", $p0, $pd) }
rep("chi(M_2) from Hurwitz numbers", (RefHur::xi0([], 2, 2, 'char'))[0], BR(-1) / 240);
rep("chi(M_3) from Hurwitz numbers", (RefHur::xi0([], 3, 2, 'char'))[0], BR(1) / 1008);
print "SELFTEST: $n checks, $bad failures\n";
