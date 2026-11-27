# status_cmp61.pl -- compare exact TR one-point values (output of status_trC1.pl / status_trhalf.pl) with the
# conjectural closed formula (6.1) of proof_M_residue.md, in the TR normalisation TR = sigma^(1-h) Xi_h (c = 1).
# Usage: perl status_cmp61.pl TAU TRFILE
use strict; use warnings; use Math::BigRat;
sub R { Math::BigRat->new(shift) }
my ($tau, $file) = (R(shift), shift);
my @B = (R(1), R("1/6"), R("-1/30"), R("1/42"), R("-1/30"), R("5/66"), R("-691/2730"), R("7/6"), R("-3617/510"));
sub fact { my $n=shift; my $r=R(1); $r*=$_ for 1..$n; $r }
sub pred { my ($h, $j) = @_; my $K = $h;
  my @lS = (R(0)); push @lS, -$B[$_]/(2*$_*fact(2*$_)) for 1..$K;
  my $e1 = -$tau; my $e2 = 1+$tau; my $sig = $tau*(1+$tau);
  my @c=(R(0)); push @c, (1+$j*$e2)*$lS[$_]*($e1**(2*$_)) + (1+$j*$e1)*$lS[$_]*($e2**(2*$_)) for 1..$K;
  my @r=(R(1)); push @r, R(0) for 1..$K;
  for my $n (1..$K){ my $s=R(0); $s += $_*$c[$_]*$r[$n-$_] for 1..$n; $r[$n]=$s/$n; }
  my $p = R(1); $p *= ($j - $_) for 1..2*$h-2;
  return $p*$r[$h]/($sig**(2*$h-1)); }
open my $f, '<', $file or die "no $file";
while (my $line = <$f>) { next unless $line =~ /g=(\d+) n=1: (.*)$/; my ($h, $rest) = ($1, $2); next if $h < 2;
  my ($n, $bad, $nz) = (0, 0, 0);
  while ($rest =~ /\((\d+)\)=(\S+)/g) { my ($j, $v) = ($1, R($2)); $n++; $nz++ unless $v->is_zero;
    my $p = pred($h, R($j)); if ($p != $v) { $bad++; print "  MISMATCH h=$h j=$j TR=$v formula=$p\n"; } }
  print "tau=$tau h=$h: compared $n values of j ($nz non-zero), mismatches with (6.1): $bad\n"; }
