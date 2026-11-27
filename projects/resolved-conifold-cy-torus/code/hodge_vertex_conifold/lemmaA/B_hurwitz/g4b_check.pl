#!/usr/bin/perl
# Hurwitz (k = 0) component of the genus-4, r = 3..5 Hodge values stored in HVC\hvc_multiprobe_g4b_out.txt (main session, exact in x),
# compared with the values from connected Hurwitz numbers + the finite formula (B2), read from runs/b2pt_4_*.out.
use strict; use warnings; use lib '.'; use HB; use Data; use Math::BigRat;
sub BR { Math::BigRat->new(shift) }
open my $fh, '<', '../../hvc_multiprobe_g4b_out.txt' or die; my @lines = <$fh>; close $fh;
my %hur; for my $f (glob("runs/b2pt_4_*.out")) { open my $g, '<', $f or next; my $l = <$g>; close $g; next unless defined $l;
  if ($l =~ /R=\(([\d,]+)\).*Hurwitz\/B2 (\S+)/) { $hur{$1} = $2 } }
my ($n, $bad) = (0, 0);
for my $l (@lines) { next unless $l =~ /^genus 4, R = \{([\d, ]+)\}, D = (\d+).*?: (0|non-zero, (.*?))\s+\[/; my ($Rs, $D, $what, $expr) = ($1, $2, $3, $4);
  my @R = split /,\s*/, $Rs; my $r = @R; my $key = join(",", @R);
  my $val;
  if ($what eq '0') { $val = BR(0) }
  else { my $ast = HB::fold(HB::parse($expr)); $val = Data::components($ast, 4, $r, {})->[0] }
  my $h = $hur{$key}; my $cmp = defined $h ? (BR($h) == $val ? "ok" : "MISMATCH") : "(no Hurwitz run)";
  $n++ if defined $h; $bad++ if $cmp eq "MISMATCH";
  printf "genus 4 R=(%s) D=%d: k=0 component of the Hodge value %s | Hurwitz/B2 %s %s\n", $key, $D, $val, $h // "-", $cmp }
print "g4b_check: $n comparisons, $bad mismatches\n";
