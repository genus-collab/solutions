#!/usr/bin/perl
# Writes the referee's job list (one line of ref_point.pl arguments per job), heavy jobs first.
# Systematic part: genus 1..3, all R with parts >= 2, r >= 3, J in a range the author did not use (author: J <= 9),
# and delta != 0 or D in the Theorem-A range.  Targeted part: genus 4, 5, 6 points not computed by the author.
use strict; use warnings;
my @jobs;
my %Jr = (1 => [10, 12], 2 => [10, 11], 3 => [10, 11]);
for my $g (1..3) { my ($lo, $hi) = @{$Jr{$g}};
  my @Rs; my $gen; $gen = sub { my ($pre, $min, $left) = @_; my $s = 0; $s += $_ for @$pre;
    push @Rs, [@$pre] if @$pre >= 3 && $s >= $lo; for my $j ($min..$left) { $gen->([@$pre, $j], $j, $left - $j) } };
  $gen->([], 2, $hi);
  for my $R (@Rs) { my $J = 0; $J += $_ for @$R; my $r = @$R; next if $r > 6; push @jobs, [$g, @$R] } }
push @jobs, map { [split ' '] } (
  "4 2 2 6", "4 2 3 5", "4 2 4 4", "4 3 3 4", "4 2 2 3 4", "4 2 3 3 3", "4 2 2 2 3 3", "4 2 2 2 2 2 3", "4 3 3 5", "4 2 2 2 6",
  "5 2 3", "5 3 4", "5 4 4", "5 2 6", "5 3 5", "5 2 2 4", "5 2 3 4", "5 2 2 2 2 2", "5 2 3 3 3", "5 2 9", "5 5 6", "5 2 2 8", "5 13",
  "6 7", "6 10", "6 12", "6 13", "6 2 2", "6 5 5", "4 2 2 2 2 2 2 2");
sub cost { my $j = shift; my ($g, @R) = @$j; my $J = 0; $J += $_ for @R; return $J + 3 * $g - 3 + 2 * @R }
@jobs = sort { cost($b) <=> cost($a) } @jobs;
print join(" ", @$_), "\n" for @jobs;
