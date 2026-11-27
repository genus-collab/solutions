#!/usr/bin/perl
# Referee: one point Xi^(0)_g(R) from the referee's Hurwitz engine (characters + E-logarithm), extracted twice
# (Lambert fit with EXTRA overdetermining orders, and Proposition D), versus the referee's (REC) and, for r = 1, Theorem B.
# Usage: perl ref_point.pl [-cj] [-x EXTRA] g j1 j2 ...
use strict; use warnings; use FindBin; use lib $FindBin::Bin; use RefHur; use RefRec; use Math::BigRat;
my $method = 'char'; my $extra = 2;
while (@ARGV && $ARGV[0] =~ /^-/) { my $o = shift; if ($o eq '-cj') { $method = 'cj' } elsif ($o eq '-x') { $extra = shift } }
my ($g, @R) = @ARGV; my $t0 = time;
my $r = @R; my $J = 0; $J += $_ for @R; my $D = $J - $r; my $delta = $D - 2 * $g + 2;
my ($p0, $ok, $pd) = RefHur::xi0(\@R, $g, $extra, $method);
my $rec = $delta == 0 ? undef : RefRec::xi($g, @R);
my $tb = $r == 1 ? RefRec::thmB_form1($g, $R[0]) : undef;
my @flags;
push @flags, $ok ? "lambert-structure-ok(+$extra)" : "LAMBERT-STRUCTURE-FAILED";
push @flags, $p0 == $pd ? "propD-agrees" : "PROPD-DISAGREES($pd)";
push @flags, defined $rec ? ($rec == $p0 ? "REC-agrees" : "REC-MISMATCH($rec)") : "REC-silent(delta=0)";
push @flags, $tb == $p0 ? "ThmB-agrees" : "THMB-MISMATCH($tb)" if defined $tb;
my $range = $D <= 2 * $g - 3 ? "in-range" : ($D == 2 * $g - 1 ? "D=2g-1" : "");
push @flags, ($p0->is_zero ? "zero-as-predicted" : "NONZERO-IN-THEOREM-A-RANGE") if $range;
printf "g=%d R=(%s) D=%d delta=%d %s: Xi0 = %s | %s [%d s, %s]\n", $g, join(",", @R), $D, $delta, $range, $p0, join(" ", @flags), time - $t0, $method;
