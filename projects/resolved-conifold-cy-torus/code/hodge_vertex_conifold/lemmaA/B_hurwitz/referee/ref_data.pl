#!/usr/bin/perl
# Referee check against the exact Hodge data (Mathematica, no Hurwitz theory, no TR), with the referee's own parser and
# negative sample points x = -2, -3, ... for the k = 0 extraction:
#  (a) ELSV normalisation: k = 0 part of mv_g<g>_l<l> versus A_g(m) from the referee's Hurwitz engine (points the author did not list
#      are included: l = 3 with large entries);
#  (b) k = 0 part of q_g<g>_r<r> at s = -1 (= Xi^(0)) versus the referee's (REC) (delta != 0), Theorem B (r = 1, all j), and the
#      referee's Hurwitz/Lambert value (selected points, including delta = 0);
#  (c) hvc_multiprobe_g4b_out.txt (genus 4, r = 3..5) k = 0 parts.
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefHur; use RefRec; use RefData; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
my $DATA = "$FindBin::Bin/../../data"; my $HVC = "$FindBin::Bin/../../..";
my ($n, $bad) = (0, 0); sub rep { my ($w, $a, $b) = @_; $n++; my $ok = (defined $a && defined $b && $a == $b); $bad++ unless $ok;
  printf "%s: %s vs %s %s\n", $w, $a // 'undef', $b // 'undef', $ok ? "ok" : "MISMATCH" }
my $part = shift // 'abc';
if ($part =~ /a/) {
  for my $spec ([1, 3, [[1,1,5], [2,3,5], [3,4,4]]], [2, 3, [[2,5,6], [1,1,8], [4,4,4]]], [3, 2, [[1,7], [6,6], [2,9]]], [3, 3, [[2,2,2], [1,2,7], [3,3,4]]],
                [4, 2, [[1,9], [5,6], [2,2]]], [4, 3, [[2,2,2], [1,1,6], [2,3,4], [6,6,6]]]) {
    my ($g, $l, $pts) = @$spec; my $ast = RefData::read_table("$DATA/mv_g${g}_l${l}.txt")->[0]{ast};
    for my $m (@$pts) { my %env; $env{"m$_"} = $m->[$_ - 1] for 1..$l; my $dv = RefData::k0($ast, $g, $l, \%env);
      my @R = grep { $_ >= 2 } @$m; my $c = grep { $_ == 1 } @$m; my $A = RefHur::A_list(\@R, $g, $c, 'char')->[$c];
      rep("(a) ELSV g=$g m=(@$m): Hodge k=0 part vs referee Hurwitz engine", $dv, $A) } } }
if ($part =~ /b/) {
  for my $spec (['q_g2_r1.txt', 2, 1, 14], ['q_g2_r2.txt', 2, 2, 13], ['q_g3_r1.txt', 3, 1, 16], ['q_g3_r2.txt', 3, 2, 13], ['q_g4_r1.txt', 4, 1, 18], ['q_g4_r2.txt', 4, 2, 13]) {
    my ($file, $g, $r, $JMAX) = @$spec; my $ast = RefData::read_table("$DATA/$file")->[0]{ast};
    my @Rs = $r == 1 ? (map { [$_] } 2..$JMAX) : (map { my $a = $_; map { [$a, $_] } $a..($JMAX - $a) } 2..int($JMAX / 2));
    for my $R (@Rs) { my %env = (s => -1); $env{"j$_"} = $R->[$_ - 1] for 1..$r; my $dv = RefData::k0($ast, $g, $r, \%env);
      my $J = 0; $J += $_ for @$R; my $delta = $J - $r - 2 * $g + 2;
      if ($delta != 0) { rep("(b) g=$g R=(@$R) delta=$delta: Hodge k=0 vs referee REC", $dv, RefRec::xi($g, @$R)) }
      else { print "(b) g=$g R=(@$R) delta=0: Hodge k=0 value $dv\n" }
      if ($r == 1) { rep("(b) g=$g j=$R->[0]: Hodge k=0 vs Theorem B (form 1)", $dv, RefRec::thmB_form1($g, $R->[0])) } } }
  # delta = 0 values and a few others by the referee's Hurwitz engine (Lambert fit)
  for my $spec (['q_g2_r1.txt', 2, [[3], [7]]], ['q_g2_r2.txt', 2, [[2,2], [2,5]]], ['q_g3_r1.txt', 3, [[5]]], ['q_g3_r2.txt', 3, [[2,4], [3,3]]], ['q_g4_r2.txt', 4, [[2,6]]]) {
    my ($file, $g, $pts) = @$spec; my $ast = RefData::read_table("$DATA/$file")->[0]{ast};
    for my $R (@$pts) { my $r = @$R; my %env = (s => -1); $env{"j$_"} = $R->[$_ - 1] for 1..$r; my $dv = RefData::k0($ast, $g, $r, \%env);
      my ($p0, $ok) = RefHur::xi0($R, $g, 2, 'char'); $n++; $bad++ unless $ok; print "(b) Lambert structure g=$g R=(@$R): ", ($ok ? "ok" : "FAILED"), "\n";
      rep("(b) g=$g R=(@$R): Hodge k=0 vs referee Hurwitz/Lambert", $dv, $p0) } } }
if ($part =~ /c/) {
  open my $fh, '<', "$HVC/hvc_multiprobe_g4b_out.txt" or die; while (my $l = <$fh>) {
    next unless $l =~ /^genus 4, R = \{([\d, ]+)\}, D = (\d+).*?: (0|non-zero, (.*?))\s+\[/; my ($Rs, $D, $what, $expr) = ($1, $2, $3, $4);
    my @R = split /,\s*/, $Rs; my $val = $what eq '0' ? BR(0) : RefData::k0(RefData::parse($expr), 4, scalar(@R), {});
    my $J = 0; $J += $_ for @R; my $delta = $J - @R - 6;
    if ($delta != 0) { rep("(c) g4b g=4 R=(@R) D=$D: Hodge k=0 vs referee REC", $val, RefRec::xi(4, @R)) }
    else { print "(c) g4b g=4 R=(@R) D=$D (delta=0): Hodge k=0 value $val\n" } } close $fh }
print "REF_DATA($part): $n checks, $bad failures\n";
