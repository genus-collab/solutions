# check_stream.pl h r : streaming check of Lemma A (M_irr(h)) from engine_h<h>_r<r>.txt without loading the polynomial:
# evaluates every component Xi^(k)_h at every multiset R of integers >= 2 with D(R) <= 2h-2 (D <= 2h-3: Lemma A points,
# D = 2h-2: control points), term by term while reading the file.
use strict; use warnings; use FindBin; use Math::BigRat;
my ($h, $r) = @ARGV;
my @pts; my $rec; $rec = sub { my ($pre, $min, $D, $k) = @_; if ($k == 0) { push @pts, [@$pre]; return }
  for (my $j = $min; ($j - 1) + $D <= 2*$h - 2; $j++) { $rec->([@$pre, $j], $j, $D + $j - 1, $k - 1); } };
$rec->([], 2, 0, $r);
open my $fh, '<', "$FindBin::Bin/engine_h${h}_r${r}.txt" or die;
my %val;   # "k|pointindex" -> BigRat
while (my $line = <$fh>) { chomp $line; next unless $line =~ /\S/;
  my ($hh, $rr, $k, $rest) = split /\s*\|\s*/, $line, 4; $rest //= '';
  my @acc = map { Math::BigRat->new(0) } @pts;
  for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $cf) = split /:/, $t; my @e = split /,/, $m;
    my $c; for my $pi (0..$#pts) { my $p = $pts[$pi]; my $mon = Math::BigInt->new(1); for my $i (0..$r-1) { $mon *= Math::BigInt->new($p->[$i]) ** $e[$i] if $e[$i]; }
      $c //= Math::BigRat->new($cf); $acc[$pi] += $c * $mon; } }
  $val{"$k|$_"} = $acc[$_] for 0..$#pts; }
close $fh;
my ($nA, $badA, $nC, $nzC) = (0, 0, 0, 0);
for my $pi (0..$#pts) { my $p = $pts[$pi]; my $D = 0; $D += $_ - 1 for @$p; my $allz = 1; my @nzk;
  for my $k (0..$h) { if (($val{"$k|$pi"} // 0) != 0) { $allz = 0; push @nzk, $k; } }
  if ($D <= 2*$h - 3) { $nA++; unless ($allz) { $badA++; print "  VIOLATION: R = {@$p}, non-zero components @nzk\n"; } }
  else { $nC++; $nzC++ unless $allz; print "  control R = {@$p} (D = $D): non-zero components {@nzk}\n"; } }
print "genus $h, r = $r: Lemma A points $nA, violations $badA; control points (D = 2h-2) $nC, non-zero $nzC\n";
