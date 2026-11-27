# compare.pl fileA fileB : compare component files (lines "h | r | k | e1,..:coef ; ...") entry by entry
use strict; use warnings; use Math::BigRat;
sub load { my $f = shift; my %d; open my $fh, '<', $f or die "$f"; while (<$fh>) { chomp; next unless /\S/;
    my ($h, $r, $k, $rest) = split /\s*\|\s*/, $_, 4; $rest //= '';
    for my $t (split /\s*;\s*/, $rest) { next unless $t =~ /:/; my ($m, $c) = split /:/, $t; $d{"$h|$r|$k|$m"} = Math::BigRat->new($c); } }
  return \%d; }
my ($A, $B) = map { load($_) } @ARGV[0,1];
my %hr; for (keys %$A, keys %$B) { my ($h, $r) = split /\|/; $hr{"$h|$r"} = 1; }
for my $x (sort keys %hr) { my ($bad, $n) = (0, 0);
  my %keys = map { $_ => 1 } grep { /^\Q$x\E\|/ } (keys %$A, keys %$B);
  for my $k (keys %keys) { $n++; my $a = $A->{$k} // Math::BigRat->new(0); my $b = $B->{$k} // Math::BigRat->new(0); $bad++ if $a != $b; }
  print "(h|r) = $x: $n coefficients compared, $bad mismatches\n"; }
