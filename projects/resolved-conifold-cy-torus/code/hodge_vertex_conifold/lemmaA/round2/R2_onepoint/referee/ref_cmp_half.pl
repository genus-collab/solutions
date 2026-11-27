#!/usr/bin/perl
# ref_cmp_half.pl TRHALF_OUT [PIDX] : compare the referee's formal Theorem 1 (ref_xy_mod.pl at tau = -1/2, no P2 needed) with the
# exact EO recursion on C_{-1/2} (status_trhalf.pl output: value = [x^mu] F_{g,n} / prod Pi_{mu_i}, odd parts only), for n = 1, 2.
# ref_xy_mod.pl prints Xi = sigma^{g-1} [x^mu]F / prod Pi, so we compare ref * sigma^{1-g} with the TR value, mod p.
use strict; use warnings;
my ($file, $PIDX) = @ARGV; $PIDX //= 0;
sub isprime { my $q = shift; for (my $d = 2; $d * $d <= $q; $d++) { return 0 if $q % $d == 0; } return 1; }
my @PR; { my $q = 2147483647; while (@PR < 4) { push @PR, $q if isprime($q); $q -= 2; } }
my $P = $PR[$PIDX];
sub mm { return ($_[0] * $_[1]) % $P; }
sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
sub pw { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub iv { return pw($_[0], $P - 2); }
sub decmod { my $s = shift; my $neg = ($s =~ s/^-//); my $r = 0; for my $ch (split //, $s) { $r = ($r * 10 + ord($ch) - 48) % $P; } return $neg ? ng($r) : $r; }
sub ratmod { my $s = shift; if ($s =~ m{^(-?\d+)/(\d+)$}) { return mm(decmod($1), iv(decmod($2))); } return decmod($s); }
my $sig = ratmod("-1/4");
open my $fh, '<', $file or die; my @lines = <$fh>; close $fh;
my ($tot, $bad) = (0, 0);
for my $line (@lines) { next unless $line =~ /^tau=-1\/2 g=(\d+) n=(\d+): (.*)$/; my ($g, $n, $rest) = ($1, $2, $3); next if $n > 2;
  my %tr; while ($rest =~ /\(([\d,]+)\)=(-?[\d\/]+)/g) { my $k = join(',', sort { $b <=> $a } split /,/, $1); $tr{$k} = ratmod($2); }
  my $jmax = 0; for my $k (keys %tr) { for (split /,/, $k) { $jmax = $_ if $_ > $jmax; } }
  my @out = `perl ref_xy_mod.pl -1/2 $g $n $jmax - $PIDX`;
  my ($c, $b) = (0, 0);
  for my $o (@out) { next unless $o =~ /mu=\(([\d,]+)\): ref (\d+)/; my ($k, $v) = ($1, $2); next unless exists $tr{$k};
    my $mine = mm($v, pw(iv($sig), $g - 1)); $c++; if ($mine != $tr{$k}) { $b++; print "  MISMATCH g=$g n=$n mu=($k): ref $mine TR $tr{$k}\n"; } }
  my $miss = scalar(keys %tr) - $c;
  print "tau=-1/2 g=$g n=$n: compared $c, mismatches $b" . ($miss ? ", not compared $miss (Pi = 0 or out of range)" : "") . "\n"; $tot += $c; $bad += $b; }
print "TOTAL tau=-1/2 (Theorem 1 formal vs exact EO, no P2): compared $tot, mismatches $bad\n";
