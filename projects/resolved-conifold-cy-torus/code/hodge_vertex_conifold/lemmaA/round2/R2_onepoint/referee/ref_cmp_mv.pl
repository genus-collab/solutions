#!/usr/bin/perl
# ref_cmp_mv.pl : compare the referee's formal Theorem 1 + P2 (ref_xy_mod.pl) with Route A's exact Marino-Vafa (Schur-side,
# Laguerre-transform) values Xi_h(R) printed in A_mv/a2_g*_r*.out (lines "h=H tau=T R=(...) D=d: Xi=V"), mod p.
use strict; use warnings;
my $P = 2147483647;
sub mm { return ($_[0] * $_[1]) % $P; } sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
sub pw { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub decmod { my $s = shift; my $neg = ($s =~ s/^-//); my $r = 0; for my $ch (split //, $s) { $r = ($r * 10 + ord($ch) - 48) % $P; } return $neg ? ng($r) : $r; }
sub ratmod { my $s = shift; if ($s =~ m{^(-?\d+)/(\d+)$}) { return mm(decmod($1), pw(decmod($2), $P - 2)); } return decmod($s); }
my $dir = "../../../A_mv";
my %want;   # "tau|h|r" -> { "R" => value }
for my $f (glob("$dir/a2_g*_r*.out")) { open my $fh, '<', $f or die; while (<$fh>) {
    next unless /^h=(\d+) tau=(\S+) R=\(([\d,]+)\) D=\d+: Xi=(-?[\d\/]+)/; my ($h, $t, $R, $v) = ($1, $2, $3, $4);
    my @R = sort { $b <=> $a } split /,/, $R; $want{"$t|$h|" . scalar(@R)}{join(',', @R)} = $v; } close $fh; }
my ($tot, $bad) = (0, 0);
for my $k (sort keys %want) { my ($t, $h, $r) = split /\|/, $k; my $jmax = 0; for my $R (keys %{$want{$k}}) { for (split /,/, $R) { $jmax = $_ if $_ > $jmax; } }
  my @out = `perl ref_xy_mod.pl $t $h $r $jmax -`; my %mine; for (@out) { $mine{$1} = $2 if /mu=\(([\d,]+)\): ref (\d+)/; }
  my ($c, $b) = (0, 0);
  for my $R (sort keys %{$want{$k}}) { $c++; my $v = ratmod($want{$k}{$R}); if (!exists $mine{$R} || $mine{$R} != $v) { $b++; print "  MISMATCH tau=$t h=$h R=($R): MV $want{$k}{$R}\n"; } }
  print "tau=$t h=$h r=$r: compared $c with Route A's MV values, mismatches $b\n"; $tot += $c; $bad += $b; }
print "TOTAL vs Route A MV (Schur side): compared $tot, mismatches $bad\n";
