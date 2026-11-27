#!/usr/bin/perl
# Referee: Theorem H tested with the independent Hurwitz-number pipeline ref_elsv.pl.
#  mode 'tables' : compare Xi^(0)_h(R) (ELSV/Burnside) with the author's k = 0 component polynomial comp_g<h>_r<r>.txt
#                  at every lattice point with parts >= 1 and D <= 2h-1 (range, boundary D = 2h-2, and D = 2h-1);
#  mode 'new'    : Theorem H at points NOT covered by the author's component tables (r >= 3, genus 5, 6), plus boundary values.
# Usage: perl ref_elsv_tests.pl tables | new5 | new6 | new3r
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use Math::BigRat; use rpoly;
require "$FindBin::Bin/ref_elsv.pl";
my $mode = shift // 'tables';
sub pts { my ($h, $r, $Dmax) = @_; my @res; my @cur; my $rec; $rec = sub { my ($i, $min) = @_; my $D = 0; $D += $_ - 1 for @cur; return if $D > $Dmax;
    if ($i == $r) { push @res, [@cur]; return; } for (my $j = $min; $j <= $Dmax + 1; $j++) { push @cur, $j; $rec->($i + 1, $j); pop @cur; } }; $rec->(0, 1); return @res; }
if ($mode eq 'tables') {
  my ($nt, $nbad) = (0, 0);
  for my $hr ([2, 1], [2, 2], [3, 1], [3, 2], [4, 1], [4, 2]) { my ($h, $r) = @$hr;
    open(my $fh, '<', "$FindBin::Bin/../comp_g${h}_r${r}.txt") or die; my $X0; while (<$fh>) { s/\r//; chomp; my @f = split /\s*\|\s*/; $X0 = parse_rp($f[3]) if $f[2] == 0; } close $fh;
    for my $p (pts($h, $r, 2 * $h - 1)) { my $D = 0; $D += $_ - 1 for @$p; my ($v) = main::xi0($h, @$p); my $w = rp_eval($X0, @$p); $nt++;
      my $ok = ($v == $w); $nbad++ unless $ok;
      printf "h=%d R=(%s) D=%d %-9s ELSV: %-22s table: %-22s %s\n", $h, join(',', @$p), $D, ($D <= 2 * $h - 3 ? 'range' : ($D == 2 * $h - 2 ? 'boundary' : 'outside')), $v->bstr, $w->bstr, ($ok ? 'agree' : 'DISAGREE'); } }
  print "SUMMARY: $nt lattice points compared, $nbad disagreements\n"; exit; }
my @list;
if ($mode eq 'new3r') { # genus 3 and 4, r = 3, 4, 5 (all parts >= 2 = irreducible points, plus their boundary D = 2h-2)
  for my $h (3, 4) { for my $r (3 .. 2 * $h - 3) { for my $p (pts($h, $r, 2 * $h - 2)) { next if grep { $_ < 2 } @$p; push @list, [$h, @$p]; } } } }
elsif ($mode eq q(new5)) { for my $j (2 .. 11) { push @list, [5, $j]; } }
elsif ($mode eq q(new5r2)) { for my $p ([2,2], [3,3], [2,7], [4,5], [3,6], [5,5], [2,8], [4,6]) { push @list, [5, @$p]; } }
elsif ($mode eq 'new5r3') { for my $r (3) { for my $p (pts(5, $r, 8)) { next if grep { $_ < 2 } @$p; push @list, [5, @$p]; } } }
elsif ($mode eq 'new6') { for my $j (2 .. 12) { push @list, [6, $j]; } }
my ($nz, $nin) = (0, 0);
for my $t (@list) { my ($h, @R) = @$t; my $D = 0; $D += $_ - 1 for @R; my $t0 = time; my ($v) = main::xi0($h, @R);
  my $cls = $D <= 2 * $h - 3 ? 'range' : ($D == 2 * $h - 2 ? 'boundary' : 'outside'); if ($cls eq 'range') { $nin++; $nz++ unless $v->is_zero; }
  printf "h=%d R=(%s) D=%d %-9s Xi^(0) = %s%s   [%d s]\n", $h, join(',', @R), $D, $cls, $v->bstr, ($cls eq 'range' ? ($v->is_zero ? '  (Theorem H: 0, OK)' : '  *** THEOREM H VIOLATED ***') : ''), time - $t0; }
print "SUMMARY ($mode): $nin points in the range, $nz non-zero\n";
