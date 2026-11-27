#!/usr/bin/perl
# Referee: genus-0 input of (REC): Xi_0(R) = (J-1)(J-2)...(J-r+3) for r >= 3 (Lambert fit of genus-0 Hurwitz numbers)
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefHur; use RefRec; use Math::BigRat;
my ($n, $bad) = (0, 0);
for my $R ([2,2,2], [2,3,5], [4,4,4], [2,2,2,2], [2,3,3,4], [2,2,2,2,2], [3,3,3,3,2], [2,2,2,2,2,2]) {
  my ($p0, $ok, $pd) = RefHur::xi0($R, 0, 3, 'char'); my $exp = RefRec::xi_genus0(@$R); $n++;
  my $good = $ok && $p0 == $exp && $pd == $exp; $bad++ unless $good;
  print "genus 0 R=(@$R): Lambert p0 $p0, PropD $pd, formula $exp ", ($good ? "ok" : "MISMATCH"), "\n" }
print "GENUS0: $n checks, $bad failures\n";
