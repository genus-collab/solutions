#!/usr/bin/perl
# compare ref_xy_mod.pl with the first-layer values of the main session's Mathematica EO recursion (hvc_tr_tau{2,3}_g3n3.txt)
use strict; use warnings;
my $P = 2147483647;
sub mm { return ($_[0] * $_[1]) % $P; } sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
sub pw { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub decmod { my $s = shift; my $neg = ($s =~ s/^-//); my $r = 0; for my $ch (split //, $s) { $r = ($r * 10 + ord($ch) - 48) % $P; } return $neg ? ng($r) : $r; }
sub ratmod { my $s = shift; if ($s =~ m{^(-?\d+)/(\d+)$}) { return mm(decmod($1), pw(decmod($2), $P - 2)); } return decmod($s); }
my ($tot, $bad) = (0, 0);
for my $t (2, 3) { open my $fh, '<', "../../../../hvc_tr_tau${t}_g3n3.txt" or die; my @L = <$fh>; close $fh;
  for my $l (@L) { next unless $l =~ /^\((\d),(\d)\):.*first layer \{(.*?)\}; TR/; my ($g, $n, $fl) = ($1, $2, $3);
    my @out = `perl ref_xy_mod.pl $t $g $n 6 -`; my %mine; for (@out) { $mine{$1} = $2 if /mu=\(([\d,]+)\): ref (\d+)/; }
    while ($fl =~ /\{([\d, ]+)\} -> (-?[\d\/]+)/g) { my $k = join(',', sort { $b <=> $a } split /,\s*/, $1); my $v = ratmod($2); $tot++;
      my $ok = (exists $mine{$k} && $mine{$k} == $v); $bad++ unless $ok; print "tau=$t g=$g mu=($k): Mathematica EO $2 -> ", ($ok ? "OK" : "MISMATCH"), "\n"; } } }
print "TOTAL first-layer values vs Mathematica EO: $tot, mismatches $bad\n";
