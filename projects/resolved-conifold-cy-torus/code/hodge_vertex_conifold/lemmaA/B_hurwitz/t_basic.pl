#!/usr/bin/perl
# Basic tests of HB / HBconn: brute-force transposition counts in S_n versus the character formula,
# genus-0 Hurwitz formula, genus-1 values, one-part ELSV.
use strict; use warnings; use lib '.'; use HB; use HBconn; use Math::BigRat; use Math::BigInt;
sub BR { Math::BigRat->new(shift) }
my $fail = 0;
# ---- brute force: count tuples (t_1..t_b) of transpositions in S_n with product = sigma^{-1}, sigma of type nu u 1^*, cycles labelled
sub perm_mul { my ($p, $q) = @_; return [ map { $p->[$q->[$_]] } 0..$#$q ] }   # (p q)(i) = p(q(i))
sub brute { my ($nu, $n, $b, $conn) = @_;
  my @tr; for my $i (0..$n-1) { for my $j ($i+1..$n-1) { my @p = 0..$n-1; @p[$i,$j] = @p[$j,$i]; push @tr, [\@p, $i, $j] } }
  # enumerate ordered tuples of disjoint labelled cycles of lengths nu: product sigma
  my @sigmas;
  my $rec; $rec = sub { my ($k, $used, $perm) = @_;
    if ($k == @$nu) { push @sigmas, [@$perm]; return }
    my $len = $nu->[$k]; my @free = grep { !$used->{$_} } 0..$n-1;
    # ordered sequences of distinct free elements of length len, modulo rotation (fix minimal first element)
    my $seq; $seq = sub { my ($s) = @_; if (@$s == $len) { my $mn = $s->[0]; for (@$s) { return if $_ < $mn }
        my @p = @$perm; for my $t (0..$len-1) { $p[$s->[$t]] = $s->[($t+1) % $len] } my %u = %$used; $u{$_} = 1 for @$s; $rec->($k+1, \%u, \@p); return }
      for my $e (@free) { next if grep { $_ == $e } @$s; $seq->([@$s, $e]) } };
    $seq->([]) };
  $rec->(0, {}, [0..$n-1]);
  my $count = 0;
  for my $sg (@sigmas) {
    my $it; $it = sub { my ($k, $prod, $edges) = @_;
      if ($k == $b) { my $tot = perm_mul($prod, $sg); for (0..$n-1) { return if $tot->[$_] != $_ }
        if ($conn) { my @par = 0..$n-1; my $find; $find = sub { my $x = shift; $x = $par[$x] while $par[$x] != $x; $x };
          for my $e (@$edges) { my ($a, $c) = @$e; $par[$find->($a)] = $find->($c) }
          for my $i (0..$n-1) { $par[$find->($i)] = $find->($sg->[$i]) }
          my $r0 = $find->(0); for (1..$n-1) { return if $find->($_) != $r0 } }
        $count++; return }
      for my $t (@tr) { $it->($k+1, perm_mul($prod, $t->[0]), [@$edges, [$t->[1], $t->[2]]]) } };
    $it->(0, [0..$n-1], []) }
  return $count }
for my $case ([[], 3, 2], [[], 3, 4], [[2], 3, 3], [[3], 3, 2], [[2], 4, 3], [[2], 4, 5], [[2,2], 4, 2], [[3], 4, 4], [[4], 4, 3], [[2], 3, 1]) {
  my ($nu, $n, $b) = @$case;
  my $bf = brute($nu, $n, $b, 0); my $T = HB::Tcount($nu, $n, [$b])->{$b};
  my $ok = $T == $bf ? "ok" : "FAIL"; $fail++ if $ok ne "ok";
  print "T(nu=[@$nu]; n=$n, b=$b): characters $T, brute force $bf  $ok\n" }
# connected counts: brute force versus exponential formula
for my $case ([[2], 3, 3], [[3], 3, 2], [[2], 4, 5], [[2,2], 4, 4], [[3], 4, 4], [[4], 4, 3], [[2], 3, 5]) {
  my ($nu, $n, $b) = @$case;
  my $bf = brute($nu, $n, $b, 1); my $C = HBconn::conn_counts($nu, $n, $b); my $v = exists $C->{$n}{$b} ? $C->{$n}{$b} : 0;
  my $ok = $v == $bf ? "ok" : "FAIL"; $fail++ if $ok ne "ok";
  print "Tc(nu=[@$nu]; n=$n, b=$b): exp. formula $v, brute force $bf  $ok\n" }
# genus 0: A_0(j,1^c) = (j+c)^(c-2)
for my $j (2..4) { my $A = HBconn::A_list([$j], 0, 5);
  for my $c (2..5) { my $exp = BR(Math::BigInt->new($j + $c) ** ($c - 2)); my $ok = $A->[$c] == $exp ? "ok" : "FAIL"; $fail++ if $ok ne "ok";
    print "A_0($j,1^$c) = $A->[$c]  expected $exp  $ok\n" } }
# genus 1: A_1(j) = (j-1)/24, A_1(j,1) = j^2/24
for my $j (2..5) { my $A = HBconn::A_list([$j], 1, 1);
  my $e0 = BR($j - 1) / 24; my $e1 = BR($j * $j) / 24;
  my $ok = ($A->[0] == $e0 && $A->[1] == $e1) ? "ok" : "FAIL"; $fail++ if $ok ne "ok";
  print "A_1($j) = $A->[0] (exp $e0), A_1($j,1) = $A->[1] (exp $e1) $ok\n" }
# genus 1, Xi^(0)_1(j) = (j-2)/24 by the finite formula
for my $j (2..6) { my $x = HBconn::xi0([$j], 1); my $e = BR($j - 2) / 24; my $ok = $x == $e ? "ok" : "FAIL"; $fail++ if $ok ne "ok";
  print "Xi0_1($j) = $x  expected $e  $ok\n" }
print "FAILURES: $fail\n";
