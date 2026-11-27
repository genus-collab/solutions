#!/usr/bin/perl
# ref_cmp_tr.pl TAU TR_OUT [NMAX] [PIDX] : compare the referee's formal Theorem 1 (ref_xy_mod.pl, no P2) with an exact EO recursion
# output whose values are [x^mu] F_{g,n} / prod Pi_{mu_i}(tau): status_trC1.pl (tau = 1, lines "g=G n=N: (mu)=v ...") or
# status_trhalf.pl (tau = -1/2, lines "tau=-1/2 g=G n=N: ..."). ref_xy_mod.pl prints sigma^{g-1} [x^mu]F / prod Pi.
use strict; use warnings;
my ($TAUS, $file, $NMAX, $PIDX) = @ARGV; $NMAX //= 9; $PIDX //= 0;
sub isprime { my $q = shift; for (my $d = 2; $d * $d <= $q; $d++) { return 0 if $q % $d == 0; } return 1; }
my @PR; { my $q = 2147483647; while (@PR < 4) { push @PR, $q if isprime($q); $q -= 2; } }
my $P = $PR[$PIDX];
sub mm { return ($_[0] * $_[1]) % $P; }
sub ng { return $_[0] == 0 ? 0 : $P - $_[0]; }
sub pw { my ($b, $e) = @_; my $r = 1; $b %= $P; while ($e > 0) { $r = ($r * $b) % $P if $e & 1; $b = ($b * $b) % $P; $e >>= 1; } return $r; }
sub iv { return pw($_[0], $P - 2); }
sub decmod { my $s = shift; my $neg = ($s =~ s/^-//); my $r = 0; for my $ch (split //, $s) { $r = ($r * 10 + ord($ch) - 48) % $P; } return $neg ? ng($r) : $r; }
sub ratmod { my $s = shift; if ($s =~ m{^(-?\d+)/(\d+)$}) { return mm(decmod($1), iv(decmod($2))); } return decmod($s); }
my $tau = ratmod($TAUS); my $sig = mm($tau, ($tau + 1) % $P);
open my $fh, '<', $file or die; my @lines = <$fh>; close $fh;
my ($tot, $bad) = (0, 0);
for my $line (@lines) { next unless $line =~ /^(?:tau=\S+ )?g=(\d+) n=(\d+): (.*)$/; my ($g, $n, $rest) = ($1, $2, $3); next if $n > $NMAX;
  my %tr; while ($rest =~ /\(([\d,]+)\)=(-?[\d\/]+)/g) { my $k = join(',', sort { $b <=> $a } split /,/, $1); $tr{$k} = ratmod($2); }
  my $jmax = 0; for my $k (keys %tr) { for (split /,/, $k) { $jmax = $_ if $_ > $jmax; } }
  my @out = `perl ref_xy_mod.pl $TAUS $g $n $jmax - $PIDX`;
  my ($c, $b) = (0, 0);
  for my $o (@out) { next unless $o =~ /mu=\(([\d,]+)\): ref (\d+)/; my ($k, $v) = ($1, $2); next unless exists $tr{$k};
    my $mine = mm($v, pw(iv($sig), $g - 1)); $c++; if ($mine != $tr{$k}) { $b++; print "  MISMATCH g=$g n=$n mu=($k): ref $mine TR $tr{$k}\n"; } }
  my $miss = scalar(keys %tr) - $c;
  print "tau=$TAUS g=$g n=$n: compared $c, mismatches $b" . ($miss ? ", not compared $miss" : "") . "\n"; $tot += $c; $bad += $b; }
print "TOTAL tau=$TAUS (Theorem 1 formal vs exact EO, no P2): compared $tot, mismatches $bad\n";
